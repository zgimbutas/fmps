function [ampole1,bmpole1]=em3mpmp3_trunc(rk,center,ampole,bmpole,nterms,center1,nterms1,radius1,rnodes,weights,nphi,ntheta)
%em3mpmp3_trunc: Apply EM multipole to multipole translation operator.
%
%  [ampole1,bmpole1]=em3mpmp3(rk,center,ampole,bmpole,nterms,...
%           center1,nterms1,radius1,rnodes,weights,nphi,ntheta);
%
%  Convert the outgoing EM multipole expansion to
%  the shifted outgoing EM multipole expansion with truncation. 
%  
%  All EM multipoles are (NTERMS+1)-by-(2*NTERMS+1) complex matrices.
%
%  Input parameters:
%
%    nterms - the number of terms in outgoing multipole expansions
%    rk - frequency parameter
%    center - real(3): sphere center location for outgoing expansion
%    ampole,bmpole - complex(ncoefs,nspheres): outgoing EM multipoles
%    center1 - real(3): sphere center location for shifted outgoing expansion
%    nterms1 - the number of terms in shifted outgoing multipole expansions
%    radius1 - real(nsource): sphere radius for the shifted outgoing expansion
%    rnodes,weights,nphi,ntheta - 
%          spherical grid, constructed via a 
%          preceding call to get_e3fgrid_trunc(nterms1,nterms1_phi)
%                 
%  Output parameters:
%
%    ampole1,bmpole1 - complex(ncoefs1): shifted outgoing EM multipoles
%
%

ncoefs1 = (nterms1+1)*(2*nterms1+1);
ampole1 = zeros(ncoefs1,1) + 1i*zeros(ncoefs1,1);
bmpole1 = zeros(ncoefs1,1) + 1i*zeros(ncoefs1,1);


mex_id_ = 'em3mpmp3trunc(i dcomplex[x], i double[], i dcomplex[], i dcomplex[], i int[x], i double[], io dcomplex[], io dcomplex[], i int[x], i double[], i double[], i double[], i int[x], i int[x])';
[ampole1, bmpole1] = fmps(mex_id_, rk, center, ampole, bmpole, nterms, center1, ampole1, bmpole1, nterms1, radius1, rnodes, weights, nphi, ntheta, 1, 1, 1, 1, 1);


