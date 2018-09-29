function [hvals,hders,emthvals,emrhvals] = emhevalrt(nterms,z)
%EMHEVALRT Spherical Bessel h functions.
%
%   [hvals,hders,emthvals,emrhvals] = emhevalrt(nterms,z);
%
%   Spherical Bessel h functions
%
%   Evaluate the h values, derivatives, tangential and radial 
%   em-scaling factors for the complex parameter z.
%
%   Input parameters:
%
%     nterms - the number of terms in em-multipole expansion
%     z - the complex argument
%
%   Output parameters:
%
%     hvals - h_n(z) values 
%     hders - h_n(z) derivatives
%     emthvals - the tangential outgoing em-scaling factors
%     emrhvals - the radial outgoing em-scaling factors
%

hvals = zeros(nterms+1,1) + 1i*zeros(nterms+1,1);
hders = zeros(nterms+1,1) + 1i*zeros(nterms+1,1);
emthvals = zeros(nterms+1,1) + 1i*zeros(nterms+1,1);
emrhvals = zeros(nterms+1,1) + 1i*zeros(nterms+1,1);

mex_id_ = 'emhevalrt(i int[x], i dcomplex[x], io dcomplex[], io dcomplex[], io dcomplex[], io dcomplex[])';
[hvals, hders, emthvals, emrhvals] = fmps(mex_id_, nterms, z, hvals, hders, emthvals, emrhvals, 1, 1);


