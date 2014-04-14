function [hvals,hders,emthvals,emrhvals] = em3hevalrt(nterms,z)
%EM3HEVALRT Spherical Bessel h functions.
%
%   [hvals,hders,emthvals,emrhvals] = em3hevalrt(nterms,z);
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

i=(0:nterms)';

hvals = besselh(i+0.5,1,z) *sqrt(pi/2/z);
hders = -(besselh(i+1.5,1,z)-i.*besselh(i+0.5,1,z)./z) *sqrt(pi/2/z);

emthvals=hvals/z+hders;
emrhvals=hvals/z.*sqrt(i.*(i+1.0d0));

% zeroth order term has no meaning for em-multipoles, set to zero
emthvals(1)=0;
emrhvals(1)=0;
