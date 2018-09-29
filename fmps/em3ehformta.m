function [aimpole,bimpole]=em3ehformta(rk, sphere_xyz, sphere_r, rnodes, weights, nphi, ntheta, evecs, hvecs, nterms)
%em3ehformta: Form incoming EM multipole from E,H values on a spherical grid.
%
%  [aimpole,bimpole]=em3ehformmp(rk, sphere_xyz, sphere_r,...
%           rnodes, weights, nphi, ntheta, evecs, hvecs, nterms);
%
%  All EM multipoles are (NTERMS+1)-by-(2*NTERMS+1) complex matrices.
%
%  Input parameters:
%
%    rk - the frequency parameter
%    sphere_xyz - real(3): sphere center location (for both E,H and multipoles)
%    sphere_r - sphere radius (for both E,H and multipoles)
%    rnodes,weights,nphi,ntheta - 
%          spherical grid, constructed via a preceding call to e3fgrid(nterms)
%    evecs,hvecs - E and H values on the spherical grid.
%    nterms - the number of terms in multipole expansions to be formed
%                 
%  Output parameters:
%
%    aimpole,bimpole - complex(ncoefs): incoming EM multipoles
%
%

ncoefs = (nterms+1)*(2*nterms+1);
aimpole = zeros(ncoefs,1) + 1i*zeros(ncoefs,1);
bimpole = zeros(ncoefs,1) + 1i*zeros(ncoefs,1);

mex_id_ = 'em3ehformta(i dcomplex[x], i double[], i double[], i dcomplex[], i dcomplex[], i double[], i double[], i int[x], i int[x], i double[x], io dcomplex[], io dcomplex[], i int[x])';
[aimpole, bimpole] = fmps(mex_id_, rk, sphere_xyz, sphere_r, evecs, hvecs, rnodes, weights, nphi, ntheta, sphere_r, aimpole, bimpole, nterms, 1, 1, 1, 1, 1);


