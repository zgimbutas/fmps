
%
%  Retrieve the reflection matrix from a file
%

filename = 'scat.freq.70';

fid=fopen(filename,'r');

s = fscanf(fid,'%s\n',1);
omega = fscanf(fid,'%g\n',1);

s = fscanf(fid,'%s\n',1);
rk = fscanf(fid,'%g',2);
rk = rk(1) + 1i*rk(2);

s = fscanf(fid,'%s\n',1);
r = fscanf(fid,'%g',1);

s = fscanf(fid,'%s\n',1);
nterms = fscanf(fid,'%d',1);

ncoefs = (nterms+1)^2;

raa = zeros(ncoefs,ncoefs)+1i*zeros(ncoefs,ncoefs);
rab = zeros(ncoefs,ncoefs)+1i*zeros(ncoefs,ncoefs);
rba = zeros(ncoefs,ncoefs)+1i*zeros(ncoefs,ncoefs);
rbb = zeros(ncoefs,ncoefs)+1i*zeros(ncoefs,ncoefs);

%
%  Loop over multipole coefficients and build the whole matrix
%
ier = 0;

for i=1:2

for k=2:ncoefs

s = fscanf(fid,'%s\n',1);
itype_mp = fscanf(fid,'%d',1);
if( i ~= itype_mp),
ier = 1;
end

s = fscanf(fid,'%s\n',1);
n = fscanf(fid,'%d',1);

s = fscanf(fid,'%s\n',1);
m = fscanf(fid,'%d',1);

s = fscanf(fid,'%s\n',2);
adata = fscanf(fid,'%g',2*ncoefs);
adata = adata(1:2:end) + 1i* adata(2:2:end);

s = fscanf(fid,'%s\n',2);
bdata = fscanf(fid,'%g',2*ncoefs);
bdata = bdata(1:2:end) + 1i* bdata(2:2:end);

if( i == 1 ),
raa(:,k) = adata;
rba(:,k) = bdata;
end

if( i == 2 ),
rab(:,k) = adata;
rbb(:,k) = bdata;
end

end

end

%
%  Zero'th term has no meaning in EM scattering, fill with zeros
%
raa(:,1) = zeros(ncoefs,1);
rab(:,1) = zeros(ncoefs,1);
rba(:,1) = zeros(ncoefs,1);
rbb(:,1) = zeros(ncoefs,1);

fclose(fid);


if( ier == 1 ),
fprintf('error in parsing data file')
end


%
%  Extract the part of reflection matrix due to dipole moments only
%
raa_dipole = raa(2:4,2:4)
rab_dipole = rab(2:4,2:4)
rba_dipole = rba(2:4,2:4) 
rbb_dipole = rbb(2:4,2:4)


%
%  Extract the part of reflection matrix due to dipole moments only
%
r_dipole= [raa(2:4,2:4) rab(2:4,2:4);
           rba(2:4,2:4) rbb(2:4,2:4)];
