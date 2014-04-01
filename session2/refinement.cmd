
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

// redefine the "r" command to adjust form factors.
r :::= {
    set vertex old_vid id;
    set edge old_eid id;

    minlen := 1.366*min(edges, length);
    // Lengths should be near min, sqrt(3)*min, and 2*min; want to catch
    // these last two cohorts.
    refine edges where length > minlen;
    u;

    foreach vertex vv where old_vid == 0 do {
        vv.ref_coord[1] :=
            avg(vv.edge ee where old_eid != 0,sum(ee.vertex where old_vid != 0, ref_coord[1]));
        vv.ref_coord[2] :=
            avg(vv.edge ee where old_eid != 0,sum(ee.vertex where old_vid != 0, ref_coord[2]));
        vv.ref_coord[3] :=
            avg(vv.edge ee where old_eid != 0,sum(ee.vertex where old_vid != 0, ref_coord[3]));
    };
    set_form_factors;
    recalc;
}
