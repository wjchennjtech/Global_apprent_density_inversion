function grav_attraction = compute_gravity_anomaly(final_density, up_boundary, low_boundary, height, nmax, nlat, nlon, R, G, Lat, Lon)
    dataro(:,1:2)=final_density(:,1:2);
    dataro(:,3) = final_density(:,3) * 1e3; %the density values were converted from g/cm³ to kg/m³
    dataup = up_boundary; datalow = low_boundary;
    data2(:,3) = dataup(:,3) * 1e3;  %the boundary units are in kilometers and need to be converted to meters
    data1(:,3) = datalow(:,3) * 1e3;  %the boundary units are in kilometers and need to be converted to meters
    datat2 = mean(data2(:,3)); 
    datat1 = mean(data1(:,3)); %calculate the average value of the boundary for subsequent use
    deltah1(:,3) = data1(:,3) - datat1;
    deltah2(:,3) = data2(:,3) - datat2;
    ratio1 = deltah1(:,3) ./ (datat1 + R);
    ratio2 = deltah2(:,3) ./ (datat2 + R);
    CCnm1=zeros(nmax+1,nmax+1);SSnm1=zeros(nmax+1,nmax+1);
    CCNM1=zeros(nmax+1,nmax+1);SSNM1=zeros(nmax+1,nmax+1);
    CCnm2=zeros(nmax+1,nmax+1);SSnm2=zeros(nmax+1,nmax+1);
    CCNM2=zeros(nmax+1,nmax+1);SSNM2=zeros(nmax+1,nmax+1);  
    CCNM = zeros(nmax+1, nmax+1); SSNM = zeros(nmax+1, nmax+1);
    data = [dataup(:,1:2), dataro(:,3)];
    [anm, bnm] = AB_matrix(data, nlat, nlon);
    [cnm, snm] = LS_SHA(anm, bnm, data, nlat, nmax);%perform spherical harmonic analysis on the density
    for n = 0:nmax
        factor = ((((R+datat2)/(R+height))^(n+3) - ((R+datat1)/(R+height))^(n+3)) * 4*pi*G*(n+1)) / (((2*n+1)*(n+3))*(1/(R+height)));
        CCNM(n+1,:) = factor * cnm(n+1,:);
        SSNM(n+1,:) = factor * snm(n+1,:);%compute the spherical harmonic coefficients of the first term in the iterative formula
    end
    for k = 1:5   % binomial expansion 
        coef1(:,1:2) = final_density(:,1:2);
        coef1(:,3) = (ratio1.^k) .* dataro(:,3);
        [anm1,bnm1] = AB_matrix(coef1,nlat,nlon);
        [cnm1,snm1] = LS_SHA(anm1,bnm1,coef1,nlat,nmax); %Perform spherical harmonic analysis through two steps
        coef2(:,1:2) = final_density(:,1:2);
        coef2(:,3) = (ratio2.^k) .* dataro(:,3);
        [anm2,bnm2] = AB_matrix(coef2,nlat,nlon);
        [cnm2,snm2] = LS_SHA(anm2,bnm2,coef2,nlat,nmax); %Perform spherical harmonic analysis through two steps  
       for n = 0:nmax
            if n+3 >= k
                C = nchoosek(n+3,k);
            else
                C = 0;
            end
            W1 = 4*pi*G*(n+1)*(R+datat1)*(((R+datat1)/(R+height))^(n+2))/((2*n+1)*(n+3));
            CCNM1(n + 1, :) = W1 * cnm1(n + 1, :) * C;
            SSNM1(n + 1, :) = W1 * snm1(n + 1, :) * C;
            W2 = 4*pi*G*(n+1)*(R+datat2)*(((R+datat2)/(R+height))^(n+2))/((2*n+1)*(n+3));
            CCNM2(n + 1, :) = W2 * cnm2(n + 1, :) * C;
            SSNM2(n + 1, :) = W2 * snm2(n + 1, :) * C;
       end
       CCnm1 = CCNM1 + CCnm1;
       SSnm1 = SSNM1 + SSnm1;
       CCnm2 = CCNM2 + CCnm2;
       SSnm2 = SSNM2 + SSnm2;%determine the spherical harmonic coefficients for the other two terms
    end   
    
    forwardCnm = CCNM - CCnm1 + CCnm2;
    forwardSnm = SSNM - SSnm1 + SSnm2;%The final spherical harmonic coefficients for forward modeling were obtained by subtracting the intermediate coefficients
    
    attraction = synthesis( forwardCnm, forwardSnm, nmax, Lat, Lon);%perform spherical harmonic synthesis using the spherical harmonic coefficients
    attraction = attraction * 1e5;
    grav_attraction = zeros(nlat*nlon, 3);
    for i = 1:nlat
        for j = 1:nlon
            grav_attraction((i-1)*nlon+j,1) = Lon(j);
            grav_attraction((i-1)*nlon+j,2) = Lat(i);
            grav_attraction((i-1)*nlon+j,3) = attraction(i,j);%output the final gravity attraction data
        end
    end
end
