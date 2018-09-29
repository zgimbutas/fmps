function [evec,hvec]=dipole3mt(rk,source,target,cmvec)
%DIPOLE3MT: Evaluate E and H fields of a magnetic dipole.
%
%   [EVEC,HVEC]=dipole3mt(RK,SOURCE,TARGET,CMVEC);
%
%   Evaluate E and H fields at the location target due
%   to the monochromatic magnetic dipole cmvec located
%   at an arbitrary source location.
%
%   Input parameters:
%
%       rk (complex *16)  - the frequency parameter
%       source (real *8 ) - the source point in R^3
%       target (real *8 ) - the target point in R^3
%       cmvec (complex *16) - the strength of the magnetic dipole
%
%   Output parameters:
%
%       evec (complex*16) - the electric field at the target
%       hvec (complex*16) - the magnetic field at the target
%

evec = zeros(3,1) + 1i*zeros(3,1);
hvec = zeros(3,1) + 1i*zeros(3,1);

mex_id_ = 'dipole3mt(i dcomplex[x], i double[x], i double[x], i dcomplex[x], io dcomplex[x], io dcomplex[x])';
[evec, hvec] = fmps(mex_id_, rk, source, target, cmvec, evec, hvec, 1, 3, 3, 3, 3, 3);




