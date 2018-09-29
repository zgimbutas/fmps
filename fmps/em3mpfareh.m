function [evecs,hvecs]=em3mpfareh(rk,center,aompole,bompole,nterms, ...
    rnodes,weights,nphi,ntheta)
%em3mpfareh: Evaluate the far field signature of the outgoing EM multipole.
%
%  [evecs,hvecs]=em3mpfareh(rk,center,aompole,bompole,nterms,...
%         sphere_xyz,sphere_r,rnodes,weights,nphi,ntheta)
%
%  All EM multipoles are (NTERMS+1)-by-(2*NTERMS+1) complex matrices.
%
%  Input parameters:
%
%    rk - the frequency parameter
%    center - center of the multipole expansion
%    aompole,bompole - complex(ncoefs): outgoing EM multipoles
%    nterms - the number of terms in multipole expansions 
%    sphere_xyz - real(3): sphere center location (for E and H fields)
%    sphere_r - sphere radius (for E and H fields)
%    rnodes,weights,nphi,ntheta - 
%          spherical grid, constructed via a preceding call to e3fgrid(nterms)
%                 
%  Output parameters:
%
%    evecs,hvecs - far field signature E and H values on the spherical grid.
%
%

nnodes=nphi*ntheta;
evecs = zeros(3,nnodes) + 1i*zeros(3,nnodes);
hvecs = zeros(3,nnodes) + 1i*zeros(3,nnodes);

mex_id_ = 'em3mpfarehfast(i dcomplex[x], i double[], i dcomplex[], i dcomplex[], i int[x], io dcomplex[], io dcomplex[], i double[], i double[], i int[x], i int[x])';
[evecs, hvecs] = fmps(mex_id_, rk, center, aompole, bompole, nterms, evecs, hvecs, rnodes, weights, nphi, ntheta, 1, 1, 1, 1);


