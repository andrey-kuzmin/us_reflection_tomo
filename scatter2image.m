function [I]=scatter2image(x_vec, y_vec, weights, x_min, x_max, y_min, y_max)
    sz_h = x_max - x_min + 1;
    sz_w = y_max - y_min + 1;
    
    I = zeros(sz_h,sz_w);
    
    x_vec = round(x_vec);
    y_vec = round(y_vec);
    
    within_region = ((x_vec >= x_min) & (x_vec <= x_max) & (y_vec >= y_min) & (y_vec <= y_max));
    
    x_vec = x_vec(within_region);
    y_vec = y_vec(within_region);
    weights_vec = weights(within_region);
    
    h_idx = y_vec - y_min + 1;
    w_idx = x_vec - x_min + 1;
    
    lin_idx = sub2ind([sz_h sz_w], h_idx, w_idx);
    
    I(lin_idx) = I(lin_idx) + weights_vec;
end