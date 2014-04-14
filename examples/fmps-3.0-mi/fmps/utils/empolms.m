function [pvec,mvec]=empolms(wavelength,rk,nspheres,aompole,bompole,nterms,center,radius,x0y0z0)
%EMPOLMS: Evaluate electric polarization and magnetization vectors.
%
%  Evaluate electric polarization and magnetization vectors 
%  for an arbitrary collection of spheres.
%
%  [pvec,mvec]=empolms(wavelength,rk,nspheres,aompole,bompole,...
%                       nterms,center,radius);
%
%  [pvec,mvec]=empolms(wavelength,rk,nspheres,aompole,bompole,...
%                       nterms,center,radius,x0y0z0);
%
%  All EM multipoles are NCOEFS = (NTERMS+1)-by-(2*NTERMS+1) complex matrices.
%
%  Input parameters:
%
%    wavelength - unused
%    rk - the frequency parameter
%    nspheres - the number of spheres
%    aompole,bompole - complex(ncoefs,nspheres): outgoing EM multipoles
%    nterms - the number of terms in multipole expansion
%    center - real(3,nspheres): sphere center locations
%    radius - real(npsheres) - sphere radii
%    x0y0z0 - real(3) - evaluate polarization and magnetization vectors 
%                       with respect to x0y0z0
%
%  Output parameter:
%
%    pvec - complex(3) - electric polarization vector
%    mvec - complex(3) - magnetic polarization (magnetization) vector
%
%

center1=[0 0 0]';

if( nargin >= 9 ),
for i=1:nspheres
  center(:,i)=center(:,i)-x0y0z0;
end
end


%%% find the size of enclosing sphere
radius1=max(sqrt(sum(center.^2,1))+radius) * 1.5;
nterms1=ceil(abs(radius1*rk)*1.2)+24;


itype1=1;
nquad1=nterms1;
nphi1=2*nquad1+1;
ntheta1=nquad1+1;

[rnodes1,weights1,nnodes1]=e3fgrid(itype1,nquad1,nphi1,ntheta1);

ncoefs1 = (nterms1+1)*(2*nterms1+1);
ampole1 = zeros(ncoefs1,1);
bmpole1 = zeros(ncoefs1,1);

for i=1:nspheres

[ampole1x,bmpole1x]=...
	em3mpmp3(rk,center(:,i),aompole(:,i),bompole(:,i),nterms,...
           center1,nterms1,radius1,rnodes1,weights1,nphi1,ntheta1);

ampole1 = ampole1 + ampole1x;
bmpole1 = bmpole1 + bmpole1x;

end

%%% finally, evaluate the polarization and magnetization vectors 
%%% due to the outgoing EM-multipole expansion

[pvec,mvec]=empol(wavelength,rk,ampole1,bmpole1,nterms1,radius1);	

