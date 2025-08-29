function [accumulated_density,gravity_rms_values] = accumulation_density_by_gradient(up_boundary, low_boundary, height, gravity_gradient, initial_density, max_iter, tol, alpha,Lon,Lat,nmax)
 
      nlat=size(Lat,2);   nlon=size(Lon,2);
    R=6378136.3;G=6.673e-11;  
    gravity_rms_values = zeros(max_iter, 1);
     accumulated_density = zeros(size(up_boundary, 1), max_iter + 2);
     cumulative_density = zeros(size(up_boundary, 1), 1);
    accumulated_density(:,1:2)=gravity_gradient(:,1:2);
    hWait = waitbar(0, 'Processing, please wait...');
    %perform gravity residual iteration
    for iter = 1:max_iter 
       fprintf("Iteration %d\n", iter);
       progress = (iter / max_iter) * 100; 
       waitbar(iter / max_iter, hWait, sprintf('%.1f%% - Iteration %d of %d', progress, iter, max_iter)); %set up a progress bar
       %perform apparent density inversion using gravity gradient data
       [pre_density, pre_rho_coeffs] = iterative_density_mapping_gradient(gravity_gradient, up_boundary, low_boundary, height, initial_density, max_iter, tol, alpha,Lon,Lat,nmax); 
        prev_gravity = gravity_gradient; 
        cumulative_density = cumulative_density + pre_density(:,3); %accumulate the density data obtained from each iterative inversion
        accumulated_density(:,iter+2) = cumulative_density; 
        updata_grav_gradient = compute_gravity_gradient(pre_density, up_boundary, low_boundary, height, nmax, nlat, nlon, R, G, Lat, Lon);%gravity forward modeling to facilitate subsequent differencing
        gravity_gradient(:,3)=prev_gravity(:,3)-updata_grav_gradient(:,3);%obtain the gravity attraction residuals
        gravityRMS=rms(gravity_gradient(:,3));
        gravity_rms_values(iter) = gravityRMS;
        fprintf("Gravity Gradient RMS Change = %.6f\n", gravityRMS);
    end       
      close(hWait)
  end