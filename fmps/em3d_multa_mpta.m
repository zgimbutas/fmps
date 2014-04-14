function [asmpole,bsmpole,A]=em3d_multa_mpta(center,radius,aompole,bompole,nterms,A)
%EM3D_MULTA_MPTAF: Direct (very-slow) EM multipole to local translation routine.
%
%  [ASMPOLE,BSMPOLE,A] = EM3D_MULTA_MPTA(CENTER,RADIUS,AOMPOLE,BOMPOLE,A);
%
%  Convert outgoing EM multipole expansions for a collection of spheres to
%  the incoming EM multipole expansions. Self interactions are NOT included.
%  
%  All EM multipoles are (NTERMS+1)-by-(2*NTERMS+1) complex matrices.
%
%  Input parameters:
%
%    center - real(3,nsource): sphere center locations
%    radius - real(nsource): sphere radii
%    aompole,bompole - complex(ncoefs,nspheres): outgoing EM multipoles
%
%    A - data structure  
%    A.nspheres - the number of spheres
%    A.nterms - the number of terms for multipole expansions
%    A.ncoefs - currently, must be set to to (NTERMS+1)*(2*NTERMS+1)
%    A.omega - angular frequency parameter
%    A.eps0 - complex: permittivity of exterior media
%    A.cmu0 - complex: permeability of exterior media
%    A.rnodes,A.weights,A.nphi,A.ntheta - 
%          spherical grid, contructed via a preceding call to e3fgrid(nterms)
%                 
%  Output parameters:
%
%    asmpole,bsmpole - complex(ncoefs,nspheres): incoming EM multipoles
%
%

rk=A.omega*sqrt(A.eps0)*sqrt(A.cmu0);

nspheres = A.nspheres;

ncoefs = (nterms+1)*(2*nterms+1);
asmpole = zeros(ncoefs,nspheres);
bsmpole = zeros(ncoefs,nspheres);

for j=1:nspheres
for i=1:nspheres

if( i == j ), continue; end

[atemp,btemp]=em3mpta3(rk,center(:,i),aompole(:,i),bompole(:,i),nterms, ...
    center(:,j),nterms,radius(j),A.rnodes,A.weights,A.nphi,A.ntheta);

asmpole(:,j) = asmpole(:,j) + atemp;
bsmpole(:,j) = bsmpole(:,j) + btemp;

end
end

