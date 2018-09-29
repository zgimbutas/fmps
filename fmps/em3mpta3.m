function [ampole1,bmpole1]=em3mpta3(rk,center,ampole,bmpole,nterms,center1,nterms1,radius1,rnodes,weights,nphi,ntheta)
%em3mpta3: Apply EM multipole to local translation operator.
%
%  [ampole1,bmpole1]=em3mpta3(rk,center,ampole,bmpole,nterms,...
%           center1,nterms1,radius1,rnodes,weights,nphi,ntheta);
%
%  Convert the outgoing EM multipole expansion to
%  the incoming EM multipole expansion. 
%  
%  All EM multipoles are (NTERMS+1)-by-(2*NTERMS+1) complex matrices.
%
%  Input parameters:
%
%    nterms - the number of terms in outgoing multipole expansions
%    rk - frequency parameter
%    center - real(3): sphere center location for outgoing expansion
%    ampole,bmpole - complex(ncoefs,nspheres): outgoing EM multipoles
%    center1 - real(3): sphere center location for incoming expansion
%    nterms1 - the number of terms in incoming multipole expansions
%    radius1 - real(nsource): sphere radius for the incoming expansion
%    rnodes,weights,nphi,ntheta - 
%          spherical grid, constructed via a preceding call to e3fgrid(nterms1)
%                 
%  Output parameters:
%
%    ampole1,bmpole1 - complex(ncoefs): incoming EM multipoles
%
%

ncoefs1 = (nterms1+1)*(2*nterms1+1);
ampole1 = zeros(ncoefs1,1) + 1i*zeros(ncoefs1,1);
bmpole1 = zeros(ncoefs1,1) + 1i*zeros(ncoefs1,1);


mex_id_ = 'em3mpta3(i dcomplex[x], i double[], i dcomplex[], i dcomplex[], i int[x], i double[], io dcomplex[], io dcomplex[], i int[x], i double[], i double[], i double[], i int[x], i int[x])';
[ampole1, bmpole1] = fmps(mex_id_, rk, center, ampole, bmpole, nterms, center1, ampole1, bmpole1, nterms1, radius1, rnodes, weights, nphi, ntheta, 1, 1, 1, 1, 1);


