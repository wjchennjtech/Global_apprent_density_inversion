function [Cnm,Snm] = LS_SHA(Anm,Bnm,data,nlat,nmax) 
   nmax1=nmax+1;
   Cnm=zeros(nmax1,nmax1);
   Snm=zeros(nmax1,nmax1);
   latmax=max(data(:,2));
   latmin=min(data(:,2));
   latinterval=(latmax-latmin)/(nlat-1);
   Lat = latmin:latinterval:latmax;
if latmin<0
   Lat = Lat+90;
end
   L = nlat;
   theta = Lat' * pi/180;
   for m = 0:L
        p  = plm(m:L, m, theta);
        ai = Anm(:, m+1);
        bi = Bnm(:, m+1);
        Cnm(m+1:L+1, m+1) = p \ ai;
        Snm(m+1:L+1, m+1) = p \ bi;  
   end
end