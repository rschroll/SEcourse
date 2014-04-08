// diagnostics.cmd

gg := { g 100; hessian_seek; }

// converge_to(max allowed relative energy step)
procedure converge_to(real relenergy) {
    recalc;
    do {
        oen := total_energy;
        gg;
        en := total_energy;
        printf "Difference: %g\n", (oen-en)/abs(en);
        printf "Target:     %g\n", relenergy;
    } while ((oen-en)/abs(en) > relenergy);
    printf "Thickness:  %g\n", thicknezz;
    printf "Gridsize:   %g = (1/%g)\n", gridsize, 1/gridsize;
    printf "Datafile:   %s\n", datafilename;
}

cc3 := {converge_to(1e-3)}
cc6 := {converge_to(1e-6)}
cc9 := {converge_to(1e-9)}
cc12 := {converge_to(1e-12)}
cc0 := {converge_to(0)}

// Output the surface - the vertices and facets, their energies, and
// the un-strained configuration
outputsurf := {
    printf "# datafilename: %s\n",datafilename;
    printf "# Surface output: vertices, faces, thickness, gridsize, length, compression \n";
    printf "%g %g %g %g %g %g 0 0\n",
        count(vertices,1), count(facets,1), thicknezz, gridsize, lngth, delta;
    printf "# Vertices: coords, bend, gbend, ref_coords\n";
    foreach vertices do
        printf "%0.15g %0.15g %0.15g %0.15g %0.15g %0.15g %0.15g %0.15g\n",
                x, y, z, bend+freeedgebend, gbend, ref_coord[1],  ref_coord[2],  ref_coord[3];
    printf "# Faces: vertices, stretch, form factors, 0\n";
    foreach facets ff do
        printf "%g %g %g %0.15g %0.15g %0.15g %0.15g 0\n",
            ff.vertex[1].id, ff.vertex[2].id, ff.vertex[3].id, stretch,
            form_factors[1],  form_factors[2],  form_factors[3];
}
