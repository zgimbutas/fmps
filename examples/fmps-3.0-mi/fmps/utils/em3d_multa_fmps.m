function [y,A]=em3d_multa_fmps(A,x)
%
%  EM multisphere solver multiplication routine
%
%  (E,H)^{scat}_{self} - R MPTA (E,H)^{scat}_{self}
%
'em3d_multa_fmps'

x0 = reshape(x,A.ncoefs,A.nspheres,2);

aimpole = x0(:,:,1);
bimpole = x0(:,:,2);

imulta_type = 2;

if( imulta_type == 1 ),
% Matlab call, very slow...
[asmpole,bsmpole]=em3d_multa_mpta(A.center,A.radius,aimpole,bimpole,A.nterms,A);
end

if( imulta_type == 2 ),
% Fortran 90, same calculation
[asmpole,bsmpole] = em3d_multa_mptaf90(A.nspheres,A.nterms,A.ncoefs,A.omega,A.eps0,A.cmu0,A.center,A.radius,aimpole,bimpole,A.rnodes,A.weights,A.nphi,A.ntheta);
end

if( imulta_type == 3 ),
% Fortran 90, same calculation with FMM
[asmpole,bsmpole] = em3d_multa_mptafmm(A.iprec,A.nspheres,A.nterms,A.ncoefs,A.omega,A.eps0,A.cmu0,A.center,A.radius,aimpole,bimpole,A.rnodes,A.weights,A.nphi,A.ntheta);
end

% Apply reflection matrices
[aompole,bompole,A]=em3d_multa_r(A.center,A.radius,asmpole,bsmpole,A.nterms,A);

atmpole = aimpole - aompole;
btmpole = bimpole - bompole;

y = reshape([atmpole btmpole], A.ncoefs*A.nspheres*2, 1);

