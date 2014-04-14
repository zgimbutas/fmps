function [evec,hvec] = em3d_mpole_targeval(nspheres,nterms,ncoefs,omega,eps0,cmu0,center,radius,aompole,bompole,ntargets,targets)
%EM3D_MPOLE_TARGEVAL: Evaluate EM multipoles at a collection of targets.
%
%  [EVEC,HVEC] = EM3D_MPOLE_TARGEVAL(NSPHERES,NTERMS,NCOEFS,...
%       OMEGA,EPS0,CMU0,CENTER,RADIUS,AOMPOLE,BOMPOLE,...
%       NTARGETS,TARGETS);
%
%  Evaluate the outgoing EM multipole expansions at a collection of targets.
%  
%  All EM multipoles are (NTERMS+1)-by-(2*NTERMS+1) complex matrices.
%
%  Input parameters:
%
%    nspheres - the number of spheres
%    nterms - the number of terms in multipole expansions
%    ncoefs - currently, must be set to to (NTERMS+1)*(2*NTERMS+1)
%    omega - angular frequency parameter
%    eps0 - complex: permittivity of exterior media
%    cmu0 - complex: permeability of exterior media
%    center - real(3,nsource): sphere center locations
%    radius - real(nsource): sphere radii
%    aompole,bompole - complex(ncoefs,nspheres): outgoing EM multipoles
%    ntargets - the number of targets
%    targets - real(3,ntargets): target locations
%                 
%  Output parameters:
%
%    evec,hvec - complex(3,ntargets): E and H fields at the targets
%
%

evec = zeros(3,ntargets) + 1i*zeros(3,ntargets);
hvec = zeros(3,ntargets) + 1i*zeros(3,ntargets);

mex_id_ = 'em3dmpoletargeval(i int[x], i int[x], i int[x], i double[x], i dcomplex[x], i dcomplex[x], i double[], i double[], i dcomplex[], i dcomplex[], i int[], i double[], io dcomplex[], io dcomplex[])';
[evec, hvec] = fmps_r2012a(mex_id_, nspheres, nterms, ncoefs, omega, eps0, cmu0, center, radius, aompole, bompole, ntargets, targets, evec, hvec, 1, 1, 1, 1, 1, 1);


