function [evec,hvec]=em3mpeval(rk,center,aompole,bompole,nterms,target)
%em3mpeval: Evaluate outgoing EM multipole expansion at a single target.
%
%  [evec,hvec]=em3mpeval(rk,center,aompole,bompole,nterms,target);
%
%  All EM multipoles are (NTERMS+1)-by-(2*NTERMS+1) complex matrices.
%
%  Input parameters:
%
%    rk - the frequency parameter
%    center - center of the multipole expansion
%    aompole,bompole - complex(ncoefs): outgoing EM multipoles
%    nterms - the number of terms in multipole expansion
%    target - real(3): target location
%                 
%  Output parameters:
%
%    evec,hvec - complex(3): E and H fields at the targets
%
%

evec = zeros(3,1) + 1i*zeros(3,1);
hvec = zeros(3,1) + 1i*zeros(3,1);

mex_id_ = 'em3mpeval(i dcomplex[x], i double[], i dcomplex[], i dcomplex[], i int[x], i double[], io dcomplex[], io dcomplex[])';
[evec, hvec] = fmps_r2012a(mex_id_, rk, center, aompole, bompole, nterms, target, evec, hvec, 1, 1);


