function [] = main_us_tomo()
        clc; close all;
        
        load input_data.mat
         
        alpha = 60;
        
        R_cm = 12.79; % tank radius, cm
        c = 1477 * 1e2; % speed of sound in water
        t_sample = 1/50e6; % sampling at 50 MHz
        spat_sample = t_sample * c; 
        
        % angular discretization, [-30,30] degrees
        angle_rng = [-30:(1/16):30]/180*pi; 
    
        traveltime_init = 6000;
        
        R_samp = round(R_cm / spat_sample);
        
        figure, imagesc(input_data),
                title('input data'), xlabel('transducer number'), ylabel('sample number');
              
        n_samp = size(input_data,1);
        n_loc = size(input_data,2);
        
        pad_sz_w = 20;
        input_data_padded = padarray(abs(input_data),[0 pad_sz_w], 'circular', 'both');
        
        %% estimate the local image orientation
        wnd_sz_mid = 7;
        [theta,~] = local_img_orientation(input_data, wnd_sz_mid, traveltime_init);
        
        [locs_x,locs_y] = transducer_locs(n_loc,R_samp);

        % reconstructed image parameters
        mult = 3.0;
             
        x_min = round(-R_samp/mult);
        x_max = round(R_samp/mult);

        y_min = round(-R_samp/mult);
        y_max = round(R_samp/mult);
        
        refl_img_x = zeros(numel(angle_rng),n_loc);
        refl_img_y = zeros(numel(angle_rng),n_loc);
        refl_img = zeros(numel(angle_rng),n_loc);
        refl_img_norm = zeros(numel(angle_rng),n_loc);
        
        sz_h = x_max - x_min + 1;
        sz_w = y_max - y_min + 1;
    
        IMG = zeros(sz_h,sz_w);
        IMG_norm = zeros(sz_h,sz_w);
        
        shapeInserter = vision.ShapeInserter('Shape','Polygons','BorderColor','White');
        
        for i_traveltime = traveltime_init:n_samp
            i_dist = i_traveltime / 2;            

            i_traveltime
            for i_angle = 1:numel(angle_rng)
                angle = angle_rng(i_angle);

                % generate the angular filter
                [filt,max_abs_h,max_abs_w] = tomo_ang_filt(angle, locs_x, locs_y, i_traveltime, n_samp, shapeInserter);

                mid_sz_h = max_abs_h + 1; 
                mid_sz_w = max_abs_w + 1; 

                c_samp =  i_traveltime;

                rng_samp = (c_samp - mid_sz_h+1):(c_samp + mid_sz_h-1);
                
                if (max(rng_samp) > n_samp)
                    break;
                end
                
                % compute the angular distribution
                corr_vec = conv2(input_data_padded(rng_samp,:),filt,'valid');

                sz_cut = pad_sz_w - mid_sz_w+1;
                corr_vec = corr_vec(sz_cut+1:end-sz_cut);

                % perform the back-projection for the current transducer location
                origin_vec = [locs_x; locs_y];
                origin_vec_norm = sqrt(origin_vec(1,:).^2 + origin_vec(2,:).^2);
                origin_vec_norm = repmat(origin_vec_norm,[2 1]);
                v_unit_vec = -origin_vec ./ origin_vec_norm;

                ROT = [cos(angle) -sin(angle); sin(angle) cos(angle)];
                v_new_matr = origin_vec + ROT*(v_unit_vec*i_dist);

                refl_img_x(i_angle,:) = v_new_matr(1,:);
                refl_img_y(i_angle,:) = v_new_matr(2,:);

                refl_img(i_angle,:) = corr_vec;
                refl_img_norm(i_angle,:) = ones(size(corr_vec)) / numel(corr_vec);
            end

            % amplify the angular distribution maxima
            [refl_img_edge] = tomo_ang_distr_softmax(theta, refl_img, input_data,angle_rng, alpha, n_loc, i_traveltime);
            
            % scatters coords to image
            IMG = IMG + scatter2image(refl_img_x(:),refl_img_y(:),refl_img_edge(:), x_min, x_max, y_min, y_max);
            IMG_norm = IMG_norm + scatter2image(refl_img_x(:),refl_img_y(:),ones(size(refl_img_norm(:))), x_min, x_max, y_min, y_max);
           
        end
        
        % plot the reflectivity field
        R_img = imresize(IMG./IMG_norm, [500 500]);
        figure, imagesc(R_img), title('reflectivity field'), axis off; 
end

