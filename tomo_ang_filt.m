function [filt,max_abs_h,max_abs_w] = tomo_ang_filt(angle, locs_x, locs_y, i_traveltime, n_samp, shapeInserter)
    n_loc = numel(locs_x);
    i_trans = round(n_loc/2);
    i_trans_mid_width = 8;
    i_dist = i_traveltime / 2;

    n_loc = numel(locs_x);
    
    ROT = [cos(angle) -sin(angle); sin(angle) cos(angle)];

    origin = [locs_x(i_trans); locs_y(i_trans)];

    v_unit = -origin / norm(origin,2);

    v_new = origin + ROT*(v_unit*i_dist);

    x_new = v_new(1);
    y_new = v_new(2);

    x_vec = repmat(x_new,[n_loc 1]);
    y_vec = repmat(y_new,[n_loc 1]);

    dist_vec = sqrt((locs_x' - x_vec).^2 + (locs_y' - y_vec).^2);
    travel_time_vec = dist_vec * 2;

    vec_h = round(travel_time_vec);
    vec_w = [1:n_loc]';
    ind_subset = (vec_h <= n_samp & vec_w >= i_trans-i_trans_mid_width & vec_w <= i_trans+i_trans_mid_width);
    vec_h = vec_h(ind_subset) - i_traveltime;
    vec_w = vec_w(ind_subset) - i_trans;
    k = numel(vec_h);

    max_abs_h = max(floor(max(abs(vec_h))),1);
    max_abs_w = floor(max(abs(vec_w)));

    I = zeros(2*max_abs_h+1,2*max_abs_w+1);

    sh_h = max_abs_h + 1;
    sh_w = max_abs_w + 1;

    poly = zeros(2*k,1);
    poly(1:2:end) = vec_w + sh_w;
    poly(2:2:end) = vec_h + sh_h;

    polygon = int32(poly);

    J = step(shapeInserter, I, polygon);
    filt = double(J);
    filt = filt / norm(filt(:),1);
end