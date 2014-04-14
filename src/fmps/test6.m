%
% Test dielectric sphere, reflection and transmission coefficients
% Shift local expansion
%

cjvec=[1;2;-3];
cmvec=[1;-2;3];

source=[1000;-2000;3000];

center=[0;0;0];
radius=20;
target=[0;-radius;0];


% material parameters

wavelength=2*pi;
omega=2*pi/wavelength;

epsE=1;
cmuE=1;
rkE=omega*sqrt(epsE)*sqrt(cmuE);

eps0=1.2;
cmu0=1;
rk0=omega*sqrt(eps0)*sqrt(cmu0);


%%%nterms=1120;
ntermsE=ceil(abs(radius*rkE)*1.2)+24
nterms0=ceil(abs(radius*rk0)*1.2)+24
nterms=max(ntermsE,nterms0)


% test local expansion

[ampole,bmpole]=em3formta(rkE,source,cjvec,cmvec,center,nterms);

'directly'
[evec1e,hvec1e]=dipole3et(rkE,source,target,cjvec);
[evec1e,hvec1e]=em3dipole3ehimp(rkE,epsE,cmuE,evec1e,hvec1e);
[evec1m,hvec1m]=dipole3mt(rkE,source,target,cmvec);
[evec1m,hvec1m]=em3dipole3ehimp(rkE,epsE,cmuE,evec1m,hvec1m);
evec1=evec1e+evec1m
hvec1=hvec1e+hvec1m

'via local'
[evec2,hvec2]=em3taeval(rkE,center,ampole,bmpole,nterms,target);
[evec2,hvec2]=em3dipole3ehimp(rkE,epsE,cmuE,evec2,hvec2)

'errors'
evec1-evec2
hvec1-hvec2

'relative errors'
norm((evec1-evec2)./evec1)
norm((hvec1-hvec2)./hvec1)


% construct reflected and transmitted field multipole expansions

r=radius;
[ra,rb]=rcoefs_die_arb21(nterms,omega,r,epsE,cmuE,eps0,cmu0);
[ta,tb]=tcoefs_die_arb21(nterms,omega,r,epsE,cmuE,eps0,cmu0);

[raa_diag,rbb_diag]=rmatr_diag(nterms,ra,rb);
  atvec = em3sphlin(nterms,ampole);
  btvec = em3sphlin(nterms,bmpole);
  aovec = atvec .* raa_diag;
  bovec = btvec .* rbb_diag;
  ampoleE = em3linsph(nterms,aovec);
  bmpoleE = em3linsph(nterms,bovec);   

[taa_diag,tbb_diag]=rmatr_diag(nterms,ta,tb);
  atvec = em3sphlin(nterms,ampole);
  btvec = em3sphlin(nterms,bmpole);
  aovec = atvec .* taa_diag;
  bovec = btvec .* tbb_diag;
  ampole0 = em3linsph(nterms,aovec);
  bmpole0 = em3linsph(nterms,bovec);   

% scattered field
[evec3,hvec3]=em3mpeval(rkE,center,ampoleE,bmpoleE,nterms,target);
[evec3,hvec3]=em3dipole3ehimp(rkE,epsE,cmuE,evec3,hvec3);

% transmitted field
[evec4,hvec4]=em3taeval(rk0,center,ampole0,bmpole0,nterms,target);
[evec4,hvec4]=em3dipole3ehimp(rk0,eps0,cmu0,evec4,hvec4);


% test jump condition, arbitrary target location on the sphere

rnorm = target-center;
rnorm = rnorm/norm(rnorm,2);

etangjump = cross(evec1+evec3,rnorm)-cross(evec4,rnorm)
htangjump = cross(hvec1+hvec3,rnorm)-cross(hvec4,rnorm)

enormjump = epsE*dot(evec1+evec3,rnorm)-eps0*dot(evec4,rnorm)
hnormjump = cmuE*dot(hvec1+hvec3,rnorm)-cmu0*dot(hvec4,rnorm)


%%%%%%%%%%%%%%%

center0=center;
radius0=radius;
nterms0=nterms;

% shift incoming expansion

center1=center+[.1;-.2;.1];
radius1=.5;
nterms1=12;

A1=get_e3fgrid(nterms1);

[ampole1,bmpole1]=em3tata3(rk0,center0,ampole0,bmpole0,nterms0,center1,nterms1,radius1,A1.rnodes,A1.weights,A1.nphi,A1.ntheta);

target1=center1+[0;0;-radius1];

'via original local'
[evec5,hvec5]=em3taeval(rk0,center0,ampole0,bmpole0,nterms0,target1);
[evec5,hvec5]=em3dipole3ehimp(rk0,eps0,cmu0,evec5,hvec5)

'via shifted local'
[evec6,hvec6]=em3taeval(rk0,center1,ampole1,bmpole1,nterms1,target1);
[evec6,hvec6]=em3dipole3ehimp(rk0,eps0,cmu0,evec6,hvec6)

'errors'
evec5-evec6
hvec5-hvec6

'relative errors'
norm((evec5-evec6)./evec5)
norm((hvec5-hvec6)./hvec5)

