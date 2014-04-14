function [ra,rb]=rcoefs_die_arb(nterms,omega,r,eps1,cmu1,eps2,cmu2)
%RCOEFS_DIE_ARB Construct reflection coefficients for a dielectric sphere.
%
%  [ra,rb]=rcoefs_die_arb(nterms,omega,r,eps1,cmu1,eps2,eps2);
%
%  Construct reflection coefficients for an arbitrary dielectric sphere.
%  Incoming wave
%
%  Input parameters:
%
%    nterms - the number of terms in multipole expansions
%    omega - angular frequency parameter
%    r - sphere radius
%    eps1 - complex: permittivity of exterior media
%    cmu1 - complex: permeability of exterior media
%    eps2 - complex: permittivity of interior media
%    cmu2 - complex: permeability of interior media
%
%  Output parameters:
%
%    ra - complex(nterms+1): TM reflection coefficients
%    rb - complex(nterms+1): TE reflection coefficients
% 


z1=omega*sqrt(eps1)*sqrt(cmu1)*r;
z2=omega*sqrt(eps2)*sqrt(cmu2)*r;


[jvals1,jders1,emtjvals1,emrjvals1] = emjevalrt(nterms,z1);
[hvals1,hders1,emthvals1,emrhvals1] = emhevalrt(nterms,z1);
[jvals2,jders2,emtjvals2,emrjvals2] = emjevalrt(nterms,z2);
[hvals2,hders2,emthvals2,emrhvals2] = emhevalrt(nterms,z2);

%
%  Zero'th term has no meaning in EM scattering
%
k = 2:(nterms+1);

%
%  Reflection coefficients: TM_21, TE_21
%
ra = -(sqrt(eps1).*sqrt(cmu2).*emtjvals1(k).*jvals2(k)-   ...
       sqrt(eps2).*sqrt(cmu1).*emtjvals2(k).*jvals1(k))./ ...
      (sqrt(eps1).*sqrt(cmu2).*emthvals1(k).*jvals2(k)-   ...
       sqrt(eps2).*sqrt(cmu1).*emtjvals2(k).*hvals1(k));

rb = -(sqrt(eps2).*sqrt(cmu1).*emtjvals1(k).*jvals2(k)-   ...
       sqrt(eps1).*sqrt(cmu2).*emtjvals2(k).*jvals1(k))./ ...
      (sqrt(eps2).*sqrt(cmu1).*emthvals1(k).*jvals2(k)-   ...
       sqrt(eps1).*sqrt(cmu2).*emtjvals2(k).*hvals1(k));

