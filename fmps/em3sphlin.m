function [cvec]=em3sphlin(nterms,ampole)
%EM3SPHLIN: Convert multipole expansion from unrolled format into linear.
%
%  [CVEC]=em3sphlin(NTERMS,AMPOLE);
%
%  AMPOLE is NCOEFS=(NTERMS+1)-by-(2*NTERMS+1) complex matrix.
%
%  Input parameters:
%
%    nterms - the number of terms in multipole expansion
%    ampole - complex(ncoefs): the EM multipole in standard unrolled format.
%
%  Output parameters:
%
%    cvec - complex((nterms+1)^2): the EM multipole in compressed linear format.
%

nvec = (nterms+1)*(nterms+1);
cvec = zeros(nvec,1) + 1i*zeros(nvec,1);

mex_id_ = 'em3sphlin(i dcomplex[], i int[x], io dcomplex[])';
[cvec] = fmps_r2012a(mex_id_, ampole, nterms, cvec, 1);


