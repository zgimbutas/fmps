function [jvals,jders,emtjvals,emrjvals] = emjevalrt(nterms,z)
%EMJEVALRT Spherical Bessel j functions.
%
%   [jvals,jders,emtjvals,emrjvals] = emjevalrt(nterms,z);
%
%   Spherical Bessel j functions
%
%   Evaluate the j values, derivatives, tangential and radial 
%   em-scaling factors for the complex parameter z.
%
%   Input parameters:
%
%     nterms - the number of terms in em-multipole expansion
%     z - the complex argument
%
%   Output parameters:
%
%     jvals - j_n(z) values 
%     jders - j_n(z) derivatives
%     emtjvals - the tangential incoming em-scaling factors
%     emrjvals - the radial incoming em-scaling factors
%

jvals = zeros(nterms+1,1) + 1i*zeros(nterms+1,1);
jders = zeros(nterms+1,1) + 1i*zeros(nterms+1,1);
emtjvals = zeros(nterms+1,1) + 1i*zeros(nterms+1,1);
emrjvals = zeros(nterms+1,1) + 1i*zeros(nterms+1,1);

mex_id_ = 'emjevalrt(i int[x], i dcomplex[x], io dcomplex[], io dcomplex[], io dcomplex[], io dcomplex[])';
[jvals, jders, emtjvals, emrjvals] = fmps(mex_id_, nterms, z, jvals, jders, emtjvals, emrjvals, 1, 1);


