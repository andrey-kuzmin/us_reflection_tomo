function [refl_img_edge] = tomo_ang_distr_softmax(theta, I_scatter, input_data,angle_rng, alpha, n_loc, i_traveltime)
    [~,angle_grid_h] = meshgrid(1:n_loc,1:numel(angle_rng));
    angle_rng_mat = repmat(angle_rng',[1 n_loc]);
    
    theta_mat = repmat(theta(i_traveltime,:),[numel(angle_rng) 1]);
    diff_angle = abs(angle_rng_mat - theta_mat);
    [~,min_idx] = min(diff_angle,[],1);

    grid_ang = abs(angle_grid_h - repmat(min_idx,[numel(angle_rng) 1]));
    
    sigma = 50.0;
    local_orient_mult = 1/2 * (1.0+exp(-grid_ang.^2/sigma^2));

    refl_img_edge = zeros(size(I_scatter,1),n_loc);
    for i_loc = 1:n_loc
        v = I_scatter(:,i_loc).*local_orient_mult(:,i_loc);

        v = exp(v*alpha) / sum(exp(v*alpha));
        max_vec = max(v);
        if (max_vec > 0.0)
          v = v / max_vec;
        end

        refl_img_edge(:,i_loc) = v * input_data(i_traveltime, i_loc);
    end
end