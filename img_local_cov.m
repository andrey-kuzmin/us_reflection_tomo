function [theta, lambda_ratio] = img_local_cov(CXX, CXY, CYY)
    C = zeros(2,2);
        
    C(1,1) = CXX;
    C(1,2) = CXY;
    C(2,1) = CXY;
    C(2,2) = CYY;
        
    [V,D] = eig(C);
        
    l1 = D(1,1);
    l2 = D(2,2);
        
    norm_fact = sqrt(l1^2 + l2^2);
        
    lambda_ratio =  (l2 / l1) * norm_fact;
        
    v1 = V(:,1);
    v2 = V(:,2);
    
    v2 = v2 * sign(v2(1));
        
    v = v2 ./ norm(v2,2);
      
    v_x = v(1);
    v_y = v(2);
    
    theta = atan2(v_y,v_x);
end