
// Set form factors so reference shape is unstretched.
set_form_factors := {
  foreach facet ff do
  { set ff.form_factors[1]
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

// redefine the "r" command to adjust form factors.
r :::= {
    set vertex old_vid id;
    set edge old_eid id;

    minlen := 1.366*min(edges, length);
    // Lengths should be near min, sqrt(3)*min, and 2*min; want to catch
    // these last two cohorts.
    refine edges where length > minlen;
    u;
    gridsize := gridsize/2;

    foreach vertex vv where old_vid == 0 do {
        vv.ref_coord[1] :=
            avg(vv.edge ee where old_eid != 0,sum(ee.vertex where old_vid != 0, ref_coord[1]));
        vv.ref_coord[2] :=
            avg(vv.edge ee where old_eid != 0,sum(ee.vertex where old_vid != 0, ref_coord[2]));
        vv.ref_coord[3] :=
            avg(vv.edge ee where old_eid != 0,sum(ee.vertex where old_vid != 0, ref_coord[3]));

        // Set correct bending energies
        if max(vv.edge, border) == 0 then {
            set vv bend;
            set vv gbend;
        } else {
            set vv freeedgebend;
        };
    };
    set_form_factors;
    set_thickness(thicknezz);
}
