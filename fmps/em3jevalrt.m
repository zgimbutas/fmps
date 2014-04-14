function [jvals,jders,emtjvals,emrjvals] = em3jevalrt(nterms,z)
%EM3JEVALRT Spherical Bessel j functions.
%
%   [jvals,jders,emtjvals,emrjvals] = em3jevalrt(nterms,z);
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

i=(0:nterms)';

if( abs(z) > 0 ),
  jvals = besselj(i+0.5,z) *sqrt(pi/2/z);
  jders = -(besselj(i+1.5,z)-i.*besselj(i+0.5,z)./z) *sqrt(pi/2/z);

  emtjvals=jvals/z+jders;
  emrjvals=jvals/z.*sqrt(i.*(i+1.0d0));

  % zeroth order term has no meaning for em-multipoles, set to zero
  emtjvals(1)=0;
  emrjvals(1)=0;
end

% z=0
if( abs(z) == 0 ),
  jvals=zeros(nterms+1,1);
  jders=zeros(nterms+1,1);
  emtjvals=zeros(nterms+1,1);
  emrjvals=zeros(nterms+1,1);
  if( nterms >= 0 ),
    jvals(1)=1;
    jders(1)=0;
    emtjvals(1)=0;
    emrjvals(1)=0;
  end
  if( nterms >= 1 ),
    jvals(2)=0;
    jders(2)=1/3.0;
    emtjvals(2)=2/3.0;
    emrjvals(2)=1/3.0*sqrt(2.0);
    jvals(3:nterms+1) = 0;
    jders(3:nterms+1) = 0;
    emtjvals(3:nterms+1) = 0;
    emrjvals(3:nterms+1) = 0;
  end
end
