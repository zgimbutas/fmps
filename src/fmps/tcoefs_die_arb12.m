function [ta,tb]=tcoefs_die_arb12(nterms,omega,r,eps1,cmu1,eps2,cmu2)
%TCOEFS_DIE_ARB12 Construct transmission coefficients for a dielectric sphere.
%
%  [ra,rb]=tcoefs_die_arb12(nterms,omega,r,eps1,cmu1,eps2,eps2);
%
%  Construct transmission coefficients for an arbitrary dielectric sphere.
%  Outgoing wave
%
%  Input parameters:
%
%    nterms - the number of terms in multipole expansions
%    omega - angular frequency parameter
%    r - sphere radius
%    eps1 - complex: permittivity of interior media
%    cmu1 - complex: permeability of interior media
%    eps2 - complex: permittivity of exterior media
%    cmu2 - complex: permeability of exterior media
%
%  Output parameters:
%
%    ta - complex(nterms+1): TM transmission coefficients
%    tb - complex(nterms+1): TE transmission coefficients
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

ima = 1i;

%
%  Transmission coefficients: TM_12, TE_12
%
ta =  -ima*eps1*(sqrt(cmu1)/sqrt(eps2))./ ...
      (sqrt(eps1).*sqrt(cmu2).*emtjvals1(k).*hvals2(k)-   ...
       sqrt(eps2).*sqrt(cmu1).*emthvals2(k).*jvals1(k));

tb =  -ima*cmu1*(sqrt(eps1)/sqrt(cmu2))./ ...
      (sqrt(eps2).*sqrt(cmu1).*emtjvals1(k).*hvals2(k)-   ...
       sqrt(eps1).*sqrt(cmu2).*emthvals2(k).*jvals1(k));


scale_ta = 1/sqrt(cmu2)*sqrt(cmu1) /r^2/omega^2 *(cmu2*eps2)/(cmu1*eps1)^2;
scale_tb = 1/sqrt(eps2)*sqrt(eps1) /r^2/omega^2 *(cmu2*eps2)/(cmu1*eps1)^2;

ta = ta*scale_ta;
tb = tb*scale_tb;

