function [ra,rb]=rcoefs_pmc(nterms,omega,r,eps1,cmu1)
%RCOEFS_PMC Construct reflection coefficients for a PMC sphere.
%
%  [ra,rb]=rcoefs_pec(nterms,omega,r,eps1,cmu1);
%
%  Construct reflection coefficients for a PMC sphere.
%  Incoming wave
%
%  Input parameters:
%
%    nterms - the number of terms in multipole expansions
%    omega - angular frequency parameter
%    r - sphere radius
%    eps1 - complex: permittivity of exterior media
%    cmu1 - complex: permeability of exterior media
%
%  Output parameters:
%
%    ra - complex(nterms+1): TM reflection coefficients
%    rb - complex(nterms+1): TE reflection coefficients
% 


z1=omega*sqrt(eps1)*sqrt(cmu1)*r;


[jvals1,jders1,emtjvals1,emrjvals1] = emjevalrt(nterms,z1);
[hvals1,hders1,emthvals1,emrhvals1] = emhevalrt(nterms,z1);

%
%  Zero'th term has no meaning in EM scattering
%
k = 2:(nterms+1);

%
%  Reflection coefficients
%

ra = -(emtjvals1(k)./emthvals1(k));

rb = -(jvals1(k)./hvals1(k));

