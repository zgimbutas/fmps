
%
%  Construct reflection coefficients for arbitrary dielectric sphere
%

%
%  Set all parameters
%

wavelength = 2*pi
omega = 2*pi/wavelength

eps1 = 1
cmu1 = 1

eps2 = 1.2
cmu2 = 1

r = 1

z1=omega*sqrt(eps1)*sqrt(cmu1)*r
z2=omega*sqrt(eps2)*sqrt(cmu2)*r

nterms = 3

[jvals1,jders1,emtjvals1,emrjvals1] = emjevalrt(nterms,z1);
[hvals1,hders1,emthvals1,emrhvals1] = emhevalrt(nterms,z1);
[jvals2,jders2,emtjvals2,emrjvals2] = emjevalrt(nterms,z2);
[hvals2,hders2,emthvals2,emrhvals2] = emhevalrt(nterms,z2);

%
%  Zero'th term has no meaning in EM scattering
%
k = 2:(nterms+1);

%
%  Reflection coefficients
%
ra = -(sqrt(eps1).*sqrt(cmu2).*emtjvals1(k).*jvals2(k)-   ...
       sqrt(eps2).*sqrt(cmu1).*emtjvals2(k).*jvals1(k))./ ...
      (sqrt(eps1).*sqrt(cmu2).*emthvals1(k).*jvals2(k)-   ...
       sqrt(eps2).*sqrt(cmu1).*emtjvals2(k).*hvals1(k))

rb = -(sqrt(eps2).*sqrt(cmu1).*emtjvals1(k).*jvals2(k)-   ...
       sqrt(eps1).*sqrt(cmu2).*emtjvals2(k).*jvals1(k))./ ...
      (sqrt(eps2).*sqrt(cmu1).*emthvals1(k).*jvals2(k)-   ...
       sqrt(eps1).*sqrt(cmu2).*emtjvals2(k).*hvals1(k))

