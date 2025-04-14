function [f]=synthesis(Cnm,Snm,nmax,Lat,Lon)
  nmnumber=0;
  for i = 1:nmax+1
     nmnumber = nmnumber + i;
  end
     cnm_snm = zeros(nmnumber, 4);
     row = 1;
  for n = 0:nmax
      for m = 0:n
        cnm_snm(row,1) = n;
        cnm_snm(row,2) = m;
        cnm_snm(row,3) = Cnm(n+1, m+1);
        cnm_snm(row,4) = Snm(n+1, m+1);
        row = row + 1;
      end
  end
   mmax=nmax;
   nlat=size(Lat,2);
   nlon=size(Lon,2);
   f=zeros(nlat,nlon);
   theta = Lat' * pi/180;
   Cnm=zeros(nmax+1,nmax+1);
   Snm=zeros(nmax+1,nmax+1);
   nmnumber=0;
   for i=1:nmax+1
       nmnumber=i+nmnumber;
   end
   for i=1: nmnumber
       Cnm(cnm_snm(i,1)+1,cnm_snm(i,2)+1)=cnm_snm(i,3);
       Snm(cnm_snm(i,1)+1,cnm_snm(i,2)+1)=cnm_snm(i,4);

   end 
   for m=0:mmax 
       n=0:nmax;
           f=plm(n,m,pi/2-theta)*Cnm(:,m+1)*cosd(m*Lon)...
            +plm(n,m,pi/2-theta)*Snm(:,m+1)*sind(m*Lon)...
            +f;
   end  
end