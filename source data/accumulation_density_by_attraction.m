function [accumulated_density,gravity_rms_values] = accumulation_density_by_attraction(up_boundary, low_boundary, height, gravity_attraction, initial_density, max_iter, tol, alpha,Lon,Lat,nmax)
    nlat=size(Lat,2);   nlon=size(Lon,2);
    R=6378136.3;G=6.673e-11;  
    gravity_rms_values = zeros(max_iter, 1);
     accumulated_density = zeros(size(up_boundary, 1), max_iter + 2);
     cumulative_density = zeros(size(up_boundary, 1), 1);
    accumulated_density(:,1:2)=gravity_attraction(:,1:2);
    hWait = waitbar(0, 'Processing, please wait...'); 
    for iter = 1:max_iter
       fprintf("Iteration %d\n", iter);
       progress = (iter / max_iter) * 100; 
       waitbar(iter / max_iter, hWait, sprintf('%.1f%% - Iteration %d of %d', progress, iter, max_iter));
       [pre_density, pre_rho_coeffs] = iterative_density_mapping_attraction(gravity_attraction, up_boundary, low_boundary, height, initial_density, max_iter, tol, alpha,Lon,Lat,nmax);
        prev_gravity = gravity_attraction; 
        cumulative_density = cumulative_density + pre_density(:,3);
        accumulated_density(:,iter+2) = cumulative_density; 
        updata_grav_attraction = compute_gravity_anomaly(pre_density, up_boundary, low_boundary, height, nmax, nlat, nlon, R, G, Lat, Lon);
        gravity_attraction(:,3)=prev_gravity(:,3)-updata_grav_attraction(:,3);
        gravityRMS=rms(gravity_attraction(:,3));
        gravity_rms_values(iter) = gravityRMS;
        fprintf("Gravity Anomaly RMS Change = %.6f\n", gravityRMS);  
    end       
    close(hWait);
  end