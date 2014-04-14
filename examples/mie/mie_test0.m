% Test Mie scattering code
% Expand incoming plane wave into the coefficients of multipole expansions.

nterms=6;

itype=1;
nquad=nterms;
nphi=2*nquad+1;
ntheta=nquad+1;

A.nterms = nterms;
A.nquad = nquad;
A.nphi = nphi;
A.ntheta = ntheta;

[A.rnodes,A.weights,A.nnodes]=e3fgrid(itype,nquad,nphi,ntheta);


nspheres=1;
center = zeros(3,nspheres);
radius = zeros(1,nspheres);
sphere_eps = zeros(1,nspheres)+1i*zeros(1,nspheres);
sphere_cmu = zeros(1,nspheres)+1i*zeros(1,nspheres);

center(1:3,1)=[0,0,0]';
radius(1)=1;
sphere_eps(1,1)=2;
sphere_cmu(1,1)=1;

omega=1;
eps0=1;
cmu0=1;
zk=omega*sqrt(eps0)*sqrt(cmu0)

%
%  STEP 4A: Define the incoming field.
%  The incoming plane wave is specified by the direction vector (kvec) and
%  the corresponding E polarization vector (epol).
%
kvec = zk*[0 0 -1];
epol = [1 0 0];

%
%  STEP 4B: get incoming E and H fields on all FMPS sphere boundaries
%

nnodes = A.nnodes;
nphi = A.nphi;
ntheta = A.ntheta;

evecs = zeros(3,nnodes,nspheres) + 1i*zeros(3,nnodes,nspheres);
hvecs = zeros(3,nnodes,nspheres) + 1i*zeros(3,nnodes,nspheres);

for j=1:nspheres

ntargets = nphi*ntheta;
targets = repmat(center(:,j),1,nnodes)+A.rnodes*radius(j);

[evecs(:,:,j),hvecs(:,:,j)] = ...
    em3d_planearb_targeval(kvec,epol,ntargets,targets);

end

%
%  STEP 4C: convert incoming E and H fields to local vector 
%           spherical harmonic expansions.
%
ncoefs = (nterms+1)*(2*nterms+1);
aimpole = zeros(ncoefs,nspheres);
bimpole = zeros(ncoefs,nspheres);

for j=1:nspheres

[aimpole(:,j),bimpole(:,j)]=...
   em3ehformta(zk, center(:,j), radius(j), ...
   A.rnodes, A.weights, A.nphi, A.ntheta,  ...
   evecs(:,:,j), hvecs(:,:,j), nterms);

end



% Reshape multipole arrays as (NTERMS+1)-by-(2*NTERMS+1) complex matrices.
amatr=reshape(aimpole(:,1),nterms+1,2*nterms+1);
bmatr=reshape(bimpole(:,1),nterms+1,2*nterms+1);

