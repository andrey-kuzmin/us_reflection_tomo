function [locs_x,locs_y] = transducer_locs(n_loc,R_samp)
    loc_delta_phi = (360 / n_loc) * (pi / 180);
    phi_loc = [1:n_loc]*loc_delta_phi;

    locs_x = R_samp * cos(phi_loc);
    locs_y = R_samp * sin(phi_loc);
end