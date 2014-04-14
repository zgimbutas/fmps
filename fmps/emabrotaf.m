function [mpout] = emabrotaf(nterms,rota,mpole)
%EMABROTF: Apply forward rotation operator to EM multipole.
%
%  MPOUT = EMABROTF(NTERMS,ROTA,MPOLE);
%
%  Apply forward rotation operator to EM multipole.
%  
%  Both MPOLE and MPOUT are (NTERMS+1)-by-(2*NTERMS+1) complex matrices.
%
%  Input parameters:
%
%    nterms - the number of terms in multipole expansion
%    rota - real(3): Euler rotation angles
%    mpole - complex(ncoefs): the EM multipole to be rotated
%                 
%  Output parameters:
%
%    mout - complex(ncoefs): the rotated EM multipole
%
%

ncoefs = (nterms+1)*(2*nterms+1);
mpout = zeros(ncoefs,1) + 1i*zeros(ncoefs,1);

mex_id_ = 'emabrotaf(i int[x], i double[], i dcomplex[], io dcomplex[])';
[mpout] = fmps_r2012a(mex_id_, nterms, rota, mpole, mpout, 1);



