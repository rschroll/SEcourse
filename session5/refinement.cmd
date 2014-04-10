// refinement.3.cmd
// Local refinement of an equilateral mesh, with freeedgebend correction.

// Set form factors from ref_coords.
set_form_factors := {
    foreach facet ff do {
        set ff.form_factors[1]
            (ff.vertex[2].ref_coord[1] - ff.vertex[1].ref_coord[1])^2
          + (ff.vertex[2].ref_coord[2] - ff.vertex[1].ref_coord[2])^2
          + (ff.vertex[2].ref_coord[3] - ff.vertex[1].ref_coord[3])^2;
        set ff.form_factors[2]
            (ff.vertex[2].ref_coord[1] - ff.vertex[1].ref_coord[1])
           *(ff.vertex[3].ref_coord[1] - ff.vertex[1].ref_coord[1])
          + (ff.vertex[2].ref_coord[2] - ff.vertex[1].ref_coord[2])
           *(ff.vertex[3].ref_coord[2] - ff.vertex[1].ref_coord[2])
          + (ff.vertex[2].ref_coord[3] - ff.vertex[1].ref_coord[3])
           *(ff.vertex[3].ref_coord[3] - ff.vertex[1].ref_coord[3]);
        set ff.form_factors[3]
            (ff.vertex[3].ref_coord[1] - ff.vertex[1].ref_coord[1])^2
          + (ff.vertex[3].ref_coord[2] - ff.vertex[1].ref_coord[2])^2
          + (ff.vertex[3].ref_coord[3] - ff.vertex[1].ref_coord[3])^2;
    }
}
set_form_factors;

procedure set_thickness(real th) {
    local pr;

    thicknezz := th;
    // Assume all facets have the same Poisson ratio
    pr := facets[1].poisson_ratio;
    bend.modulus := thicknezz**2 / 6 / (1 - pr**2);
    gbend.modulus := -bend.modulus * (1 - pr)/2;

    //edgewidth := 1/3;  // Right triangles along edge
    edgewidth := sqrt(3)/4;  // Equilateral triangles along edge
    freeedgebend.modulus := bend.modulus/4 * (1 - pr**2) * edgewidth * gridsize;
    recalc;
}

procedure refinemarked(){
    local newmark;
    local hypid;
    local baseid;
    local heightid;
    local otherhypid;
    local len;
    local onedge;

    set edge rheight 0;

    do {
        newmark := 0;
        foreach facet ff where sum(ff.edges, divide) do {
            // If it is an equilateral triangle
            if abs(ff.form_factors[1] - ff.form_factors[3]) < 0.01*ff.form_factors[1] then {
                // If 2/3 of the sides are to be divided, divide all
                if sum(ff.edges, divide) == 2 then {
                    set ff.edges divide 1;
                    newmark := 1;
                };
            } else {
                // The longest edge is the hypotenuse of the right triangle
                len := 0;
                foreach ff.edges fe do {
                    if fe.length > len then {
                        hypid := fe.id;
                        len := fe.length;
                    };
                };
                // The next longest is the height
                len := 0;
                foreach ff.edge fe where id != hypid do {
                    if fe.length > len then {
                        heightid := fe.id;
                        fe.rheight := 1;
                        len := fe.length;
                    };
                };
                // The last one is the base
                foreach ff.edge fe where id != hypid and id != heightid do {
                    baseid := fe.id;
                    fe.rheight := 0;
                };

                // Make sure the hypotenuse will be divided
                if not edges[hypid].divide then {
                    set edges[hypid] divide 1;
                    newmark := 1;
                };

                // If the height is on the border of the sheet, divide it
                if edges[heightid].border then {
                    if not edges[heightid].divide then {
                        set edges[heightid] divide 1;
                        newmark := 1;
                    };

                // Otherwise, don't divide it.  Instead divide the hypotenuse
                // of the pair triangle
                } else {
                    set edges[heightid] divide 0;

                    // Find the longest side of the triangle sharing the height
                    // edge - that is the hypotenuse of the pair triangle
                    foreach edges[heightid].facets ef where id != ff.id do {
                        len := 0;
                        foreach ef.edges efe do {
                            if efe.length > len then {
                                otherhypid := efe.id;
                                len := efe.length;
                            };
                        };
                    };
                    if not edges[otherhypid].divide then {
                        set edges[otherhypid] divide 1;
                        newmark := 1;
                    };
                };

                // If two bases are touching, or base is on perimeter, don't divide it
                if edges[baseid].divide then {
                    onedge := 1;
                    // Get other face
                    foreach edges[baseid].facets ef where id != ff.id do {
                        onedge := 0;
                        // Is it a right triangle?
                        if abs(ef.form_factors[1] - ef.form_factors[3]) > 0.01*ef.form_factors[1] then
                            // Is baseid also a base here?
                            if edges[baseid].length == min(ef.edges, length) then
                                set edges[baseid] divide 0;
                    };
                    if onedge then
                        set edges[baseid] divide 0;
                };
            };
        };
    } while newmark;

    set vertex old_vid id;
    set edge old_eid id;

    refine edges where divide;
    set edge divide 0;
    //gridsize := gridsize/2;

    foreach vertex vv where old_vid == 0 do {
        vv.ref_coord[1] :=
            avg(vv.edge ee where old_eid != 0,sum(ee.vertex where old_vid != 0,
                    ref_coord[1]));
        vv.ref_coord[2] :=
            avg(vv.edge ee where old_eid != 0,sum(ee.vertex where old_vid != 0,
                    ref_coord[2]));
        vv.ref_coord[3] :=
            avg(vv.edge ee where old_eid != 0,sum(ee.vertex where old_vid != 0,
                    ref_coord[3]));
        // Set bending energies for those not on an edge:
        if max(vv.edge, border) == 0 then {
            set vv bend;
            set vv gbend;
        }
        else
            set vv freeedgebend;
    };
    // Only swap edges that were newly added or were heights of right triangles
    do {
        flush_counts;
        equiangulate edge where rheight or old_eid==0;
    } while equi_count;

    set_form_factors;
    // Edge correction needs to be updated:
    set_thickness(thicknezz);
}

r :::= {
    set edge divide 1 where length > 0.683*gridsize;
    refinemarked();
    gridsize := gridsize/2;
}

function real ssqrt(real val) {
    if val < 0 then
        return 0;
    return sqrt(val);
}

procedure mark_curv(real curvgrid) {
    // Mark edges to be refined where curvature * local gridsize > curvgrid
    local pr;
    pr := facets[1].poisson_ratio;

    foreach vertex vv do
        if vv.bend + 1/(1-pr) * vv.gbend + sqrt(vv.bend) * ssqrt(vv.bend + 2/(1-pr) * vv.gbend)
                > bend.modulus * curvgrid**2 / 4 then
            set vv.edges divide 1;
}

procedure mark_stretch(real frac) {
    // Mark edges where jump in stretch density > frac * max energy density
    local maxen;
    maxen := max(facets, stretch/area);
    // For some reason, max(ee.facets, stretch) will report erroneous results,
    // while max(ee.facets, facets[id].stretch) works properly....
    set edge ee divide 1 where max(ee.facets, facets[id].stretch/area)
        - min(ee.facets, facets[id].stretch/area) > maxen * frac;
}

function integer refine_cs(real curvgrid, real frac) {
    mark_curv(curvgrid);
    mark_stretch(frac);
    set edge divide 0 where length < thicknezz;
    if sum(edges, divide) > 0 then {
        refinemarked();
        return 1;
    };
    return 0;
}

