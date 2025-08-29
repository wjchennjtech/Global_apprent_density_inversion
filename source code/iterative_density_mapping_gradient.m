function [final_density, final_rho_coeffs] = iterative_density_mapping_gradient(gravity_gradient,up_boundary,low_boundary,height,initial_density,max_iter,tol,alpha,Lon,Lat,nmax)
   
    nmax1=nmax;    
    nlon=size(Lon,2);nlat=size(Lat,2);R=6378136.3;G=6.673e-11; 
    data2(:,3) = up_boundary(:,3) * 1e3;  %the boundary units are in kilometers and need to be converted to meters
    data1(:,3) = low_boundary(:,3) * 1e3;  %the boundary units are in kilometers and need to be converted to meters
    datat2 = mean(data2(:,3)); %calculate the average value of the boundary for subsequent use
    datat1 = mean(data1(:,3));%calculate the average value of the boundary for subsequent use
    density_rms_values = zeros(max_iter, 1); 
    prev_density=zeros(length(data2(:,3)),3);
    prev_Cnm = zeros(nmax+1, nmax+1);
    prev_Snm = zeros(nmax+1, nmax+1);


  for iter = 1:max_iter %main iterative loop
    if iter==1
      rho(:,1:2)= up_boundary(:,1:2);
      rho(:,3) = initial_density;
      prev_density(:,3)=0;
    end
    rho(:,3) = rho(:,3) * 1000;%the density values were converted from g/cm³ to kg/m³
    
 % Step 1: Solve for the spherical harmonic coefficients of the first term (the fixed term).
    gravity_gradient(:,3) = gravity_gradient(:,3) / 1e9;  
    [Anm,Bnm] = AB_matrix(gravity_gradient,nlat,nlon);
    [hcnm,hsnm] = LS_SHA(Anm,Bnm,gravity_gradient,nlat,nmax);%perform spherical harmonic analysis on the gravity gradient
    CNM = zeros(nmax+1, size(hcnm,2));
    SNM = zeros(nmax+1, size(hsnm,2));
    for n = 0:nmax
        factor = ((2*n+1)*(n+3)) / ( (((R+datat2)/(R+height))^(n+3) - ((R+datat1)/(R+height))^(n+3) ) * 4*pi*G*(n+1)*(n+2) );
        CNM(n+1,:) = factor * hcnm(n+1,:);
        SNM(n+1,:) = factor * hsnm(n+1,:);%obtain the spherical harmonic coefficients of the fixed term
    end
% Step 2: Solve for the spherical harmonic coefficients of the other two iterative approximation terms.
     deltah1(:,3) = data1(:,3) - datat1;
     deltah2(:,3) = data2(:,3) - datat2;
     ratio1 = deltah1(:,3) ./ (datat1 + R);
     ratio2 = deltah2(:,3) ./ (datat2 + R);
     Cnm1=zeros(nmax+1,nmax+1);
     Snm1=zeros(nmax+1,nmax+1);
     CNM1=zeros(nmax+1,nmax+1);
     SNM1=zeros(nmax+1,nmax+1);
     Cnm2=zeros(nmax+1,nmax+1);
     Snm2=zeros(nmax+1,nmax+1);
     CNM2=zeros(nmax+1,nmax+1);
     SNM2=zeros(nmax+1,nmax+1);  

     for k = 1:5  % binomial expansion
        field1(:,1:2) = rho(:,1:2);
        field1(:,3) = (ratio1.^k) .* rho(:,3);
        [Anm1,Bnm1] = AB_matrix(field1,nlat,nlon);
        [cnm1,snm1] = LS_SHA(Anm1,Bnm1,field1,nlat,nmax);

        field2(:,1:2) = rho(:,1:2);
        field2(:,3) = (ratio2.^k) .* rho(:,3);
        [Anm2,Bnm2] = AB_matrix(field2,nlat,nlon);
        [cnm2,snm2] = LS_SHA(Anm2,Bnm2,field2,nlat,nmax);   

       for n = 0:nmax
            if n+3 >= k
                comb = nchoosek(n+3,k);
            else
                comb = 0;
            end
            F1 = 1 / (((R + datat2) / (R + datat1))^(n + 3) - 1);
            F2 = 1 / (1 - ((R + datat1) / (R + datat2))^(n + 3));
            CNM1(n + 1, :) = F1 * cnm1(n + 1, :) * comb;
            SNM1(n + 1, :) = F1 * snm1(n + 1, :) * comb;
            CNM2(n + 1, :) = F2 * cnm2(n + 1, :) * comb;
            SNM2(n + 1, :) = F2 * snm2(n + 1, :) * comb;
       end
       Cnm1 = CNM1 + Cnm1;
       Snm1 = SNM1 + Snm1;
       Cnm2 = CNM2 + Cnm2;
       Snm2 = SNM2 + Snm2;  %Solve for the spherical harmonic coefficients of the other two iterative approximation terms.
     end   

% In the first iteration, there are no previous spherical harmonic coefficients, so multiplication by alpha is not needed.
% From the second iteration onward, the spherical harmonic coefficients are a combination of (1 - alpha) times the previous coefficients and alpha times the current ones.
     if iter == 1
         lastCnm = CNM + Cnm1 - Cnm2;
         lastSnm = SNM + Snm1 - Snm2;
     else
         lastCnm = alpha * (CNM + Cnm1 - Cnm2) + (1-alpha) * prev_Cnm;
         lastSnm = alpha * (SNM + Snm1 - Snm2) + (1-alpha) * prev_Snm;
     end  

     prev_Cnm = lastCnm;
     prev_Snm = lastSnm; 
   
     

     [Iteration_of_density]=synthesis(lastCnm,lastSnm,nmax1,Lat,Lon);%perform spherical harmonic synthesis using the final obtained coefficients
     Iteration_of_density=Iteration_of_density/1e3;% convert to grams per cubic centimeter.
    
     fresult=zeros(nlat*nlon,3);    
     for i = 1:nlat
        for j = 1:nlon
            fresult((i-1)*nlon+j,1) = Lon(j);
            fresult((i-1)*nlon+j,2) = Lat(i);
            fresult((i-1)*nlon+j,3) = Iteration_of_density(i,j);%output the final density data
        end
     end
       


     density_rms_change = rms(fresult(:,3) - prev_density(:,3));
     density_rms_values(iter) = density_rms_change;
     fprintf("Density of RMS variation = %.6f\n", density_rms_change);

    nmnumber=0;
    for i = 1:nmax+1
        nmnumber = nmnumber + i;
    end
    rho_coeffs_final = zeros(nmnumber, 4);
    row = 1;
    for n = 0:nmax
        for m = 0:n
          rho_coeffs_final(row,1) = n;
          rho_coeffs_final(row,2) = m;
          rho_coeffs_final(row,3) = lastCnm(n+1, m+1);
          rho_coeffs_final(row,4) = lastSnm(n+1, m+1);%save the spherical harmonic coefficients for later use
          row = row + 1;
        end
    end
    

     final_rho_coeffs=rho_coeffs_final;
     final_density = fresult;


    if density_rms_change < tol
        fprintf("Density RMS change %.6e is smaller than the tolerance %.6e, iteration terminated.\n", density_rms_change, tol);
        break;
    end

    prev_density(:,3) = fresult(:,3);%the final output was used as the initial input for the next iteration
    rho(:,3) = fresult(:,3);%the final output was used as the initial input for the next iteration
  end

end