function [ampole,bmpole]=em3formta(rk,source,cjvec,cmvec,center,nterms)
%em3formta: Form the incoming EM multipole expansion.
%
%  [ampole,bmpole]=em3formta(rk,source,cjvec,cmvec,center,nterms);
%
%  All EM multipoles are (NTERMS+1)-by-(2*NTERMS+1) complex matrices.
%
%  Input parameters:
%
%    rk - the frequency parameter
%    center - center of the multipole expansion
%    source - real(3,nsource): the source locations in R^3
%    cjvec - complex(3,nsource): the strengths of the electric dipoles 
%    cmvec - complex(3,nsource): the strengths of the magnetic dipoles  
%    nterms - the number of terms in multipole expansion
%                 
%  Output parameters:
%
%    ampole,bmpole - complex(ncoefs): incoming EM multipoles
%
%

npts=size(source,2);
ncoefs = (nterms+1)*(2*nterms+1);
ampole = zeros(ncoefs,1) + 1i*zeros(ncoefs,1);
bmpole = zeros(ncoefs,1) + 1i*zeros(ncoefs,1);

mex_id_ = 'em3formta(i dcomplex[x], i double[], i dcomplex[], i dcomplex[], i int[x], i double[], io dcomplex[], io dcomplex[], i int[x])';
[ampole, bmpole] = fmps(mex_id_, rk, source, cjvec, cmvec, npts, center, ampole, bmpole, nterms, 1, 1, 1);



