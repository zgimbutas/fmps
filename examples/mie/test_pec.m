
%
%  Construct reflection coefficients for PEC sphere
%

%
%  Set all parameters
%

omega = 1

eps1 = 1
cmu1 = 1

r = 1

z1=omega*sqrt(eps1)*sqrt(cmu1)*r

nterms = 3

[jvals1,jders1,emtjvals1,emrjvals1] = emjevalrt(nterms,z1);
[hvals1,hders1,emthvals1,emrhvals1] = emhevalrt(nterms,z1);

%
%  Zero'th term has no meaning in EM scattering
%
k = 2:(nterms+1);

%
%  Reflection coefficients
%

ra = -(jvals1(k)./hvals1(k))

rb = -(emtjvals1(k)./emthvals1(k))

