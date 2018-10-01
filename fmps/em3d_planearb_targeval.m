function [evec,hvec]=em3d_planearb_targeval(rkvec,epol,ntargets,targets)
%EM3D_PLANEARB_TARGEVAL: Evaluate an arbitrary plane wave at a collection of targets.
%
%  [EVEC,HVEC] = EM3D_PLANEARB_TARGEVAL(RKVEK,EPOL,NTARGETS,TARGETS);
%
%  Evaluate an arbitrary oriented plane wave at a collection of targets.
%  
%  Input parameters:
%
%    rkvec - the vector Helmholz parameter
%    epol - the polarization vector
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

mex_id_ = 'emplanearbtargeval(i dcomplex[], i dcomplex[], i int[x], i double[], io dcomplex[], io dcomplex[])';
[evec, hvec] = fmps_r2012a(mex_id_, rkvec, epol, ntargets, targets, evec, hvec, 1);


