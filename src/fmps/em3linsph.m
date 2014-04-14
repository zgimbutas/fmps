function [ampole]=em3linsph(nterms,cvec)
%EM3LINSPH: Convert multipole expansion from linear format into unrolled.
%
%  [AMPOLE]=em3linsph(NTERMS,CVEC);
%
%  AMPOLE is NCOEFS=(NTERMS+1)-by-(2*NTERMS+1) complex matrix.
%
%  Input parameters:
%
%    nterms - the number of terms in multipole expansion
%    cvec - complex((nterms+1)^2): the EM multipole in compressed linear format.
%
%  Output parameters:
%
%    ampole - complex(ncoefs): the EM multipole in standard unrolled format.
%

ncoefs = (nterms+1)*(2*nterms+1);
ampole = zeros(ncoefs,1) + 1i*zeros(ncoefs,1);

mex_id_ = 'em3linsph(io dcomplex[], i int[x], i dcomplex[])';
[ampole] = fmps_r2012a(mex_id_, ampole, nterms, cvec, 1);




