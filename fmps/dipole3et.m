function [evec,hvec]=dipole3et(rk,source,target,cjvec)
%DIPOLE3ET: Evaluate E and H fields of an electric dipole.
%
%   [EVEC,HVEC]=dipole3et(RK,SOURCE,TARGET,CJVEC);
%
%   Evaluate E and H fields at the location target due
%   to the monochromatic electric dipole cjvec located 
%   at an arbitrary source location.
%
%   Input parameters:
%
%       rk (complex *16)  - the frequency parameter
%       source (real *8 ) - the source point in R^3
%       target (real *8 ) - the target point in R^3
%       cjvec (complex *16) - the strength of the electric dipole   
%
%   Output parameters:
%
%       evec (complex*16) - the electric field at the target
%       hvec (complex*16) - the magnetic field at the target
%

evec = zeros(3,1) + 1i*zeros(3,1);
hvec = zeros(3,1) + 1i*zeros(3,1);

mex_id_ = 'dipole3et(i dcomplex[x], i double[x], i double[x], i dcomplex[x], io dcomplex[x], io dcomplex[x])';
[evec, hvec] = fmps_r2012a(mex_id_, rk, source, target, cjvec, evec, hvec, 1, 3, 3, 3, 3, 3);


