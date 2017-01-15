function [theta_vec, lambda_ratio_vec] = local_img_orientation(img, wnd_sz_mid, traveltime_init)
    [h,w] = size(img);
    
    patch_sz = [2*wnd_sz_mid+1 2*wnd_sz_mid+1];
    
    theta_vec = zeros(h,w);
    v_x_vec = zeros(h,w);
    v_y_vec = zeros(h,w);
    lambda_ratio_vec = zeros(h,w);    
    
    wnd_sz = 2 * wnd_sz_mid + 1;
    
    n = wnd_sz * wnd_sz;
    
    [h_matr,w_matr] = ndgrid(1:wnd_sz,1:wnd_sz);
       
    r_x = w_matr - (wnd_sz_mid + 1);
    r_y = -(h_matr - (wnd_sz_mid + 1));
    
    r_x_rep = repmat(r_x(:),[1 w]);
    r_y_rep = repmat(r_y(:),[1 w]);
       
    r = sqrt(r_x.^2 + r_y.^2);
    alpha = exp(-(r/4).^2);
    
    alpha_rep = repmat(alpha(:),[1 w]);
    
    img_padded = padarray(img,[0 wnd_sz_mid],'symmetric','both');
    
    for i_h = traveltime_init:h-wnd_sz_mid
        i_h 
        
        strip = img_padded(i_h-wnd_sz_mid:i_h+wnd_sz_mid,:);
        
        patches = im2col(strip,patch_sz,'sliding');
  
        v_x_matr = patches .* r_x_rep .* alpha_rep;
        v_y_matr = patches .* r_y_rep .* alpha_rep;
        
        v_x_matr_c = bsxfun(@minus,v_x_matr,sum(v_x_matr,1)/n);
        v_y_matr_c = bsxfun(@minus,v_y_matr,sum(v_y_matr,1)/n);
        
        CXX = sum(v_x_matr_c .* v_x_matr_c,1) / n;
        CYY = sum(v_y_matr_c .* v_y_matr_c,1) / n;
        CXY = sum(v_x_matr_c .* v_y_matr_c,1) / n;
        
        [theta, lambda_ratio] = arrayfun(@img_local_cov,CXX,CXY,CYY);
        
        theta_vec(i_h,:) = theta;
        lambda_ratio_vec(i_h,:) = lambda_ratio;
    end
end



