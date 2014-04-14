%
% Test dielectric sphere, reflection and transmission coefficients
% Shift local expansion, reflect from gold sphere
% Shift multipole, reflect from outer sphere
%
% All dimensions are in nm
%

cjvec=[1;2;-3];
cmvec=[1;-2;3];

source=[10000;-20000;30000];

center=[0;0;0];
radius=2400;
target=[0;-radius;0];


% material parameters

%%%wavelength=2*pi;
wavelength=620;
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
[raE0,rbE0]=rcoefs_die_arb21(nterms,omega,r,epsE,cmuE,eps0,cmu0);
[taE0,tbE0]=tcoefs_die_arb21(nterms,omega,r,epsE,cmuE,eps0,cmu0);

[raa_diag,rbb_diag]=rmatr_diag(nterms,raE0,rbE0);
  atvec = em3sphlin(nterms,ampole);
  btvec = em3sphlin(nterms,bmpole);
  aovec = atvec .* raa_diag;
  bovec = btvec .* rbb_diag;
  ampoleE = em3linsph(nterms,aovec);
  bmpoleE = em3linsph(nterms,bovec);   

[taa_diag,tbb_diag]=rmatr_diag(nterms,taE0,tbE0);
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
radius1=50;
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


% construct reflected and transmitted field multipole expansions

% JCh gold at 620nm
eps1=(0.2039+3.3056i)^2;
cmu1=1;
rk1=omega*sqrt(eps1)*sqrt(cmu1);

[ra,rb]=rcoefs_die_arb21(nterms1,omega,radius1,eps0,cmu0,eps1,cmu1);
[ta,tb]=tcoefs_die_arb21(nterms1,omega,radius1,eps0,cmu0,eps1,cmu1);

[raa_diag,rbb_diag]=rmatr_diag(nterms1,ra,rb);
  atvec = em3sphlin(nterms1,ampole1);
  btvec = em3sphlin(nterms1,bmpole1);
  aovec = atvec .* raa_diag;
  bovec = btvec .* rbb_diag;
  aompole1 = em3linsph(nterms1,aovec);
  bompole1 = em3linsph(nterms1,bovec);   

'via mpole'
[evec7,hvec7]=em3mpeval(rk0,center1,aompole1,bompole1,nterms1,target);
[evec7,hvec7]=em3dipole3ehimp(rk0,eps0,cmu0,evec7,hvec7)


% shift outgoing expansion

A0=get_e3fgrid(nterms0);

[aompole0,bompole0]=em3mpmp3(rk0,center1,aompole1,bompole1,nterms1,center0,nterms0,radius0,A0.rnodes,A0.weights,A0.nphi,A0.ntheta);

'via original multipole'
[evec5,hvec5]=em3mpeval(rk0,center1,aompole1,bompole1,nterms1,target);
[evec5,hvec5]=em3dipole3ehimp(rk0,eps0,cmu0,evec5,hvec5)

'via shifted multipole'
[evec6,hvec6]=em3mpeval(rk0,center0,aompole0,bompole0,nterms0,target);
[evec6,hvec6]=em3dipole3ehimp(rk0,eps0,cmu0,evec6,hvec6)

'errors'
evec5-evec6
hvec5-hvec6

'relative errors'
norm((evec5-evec6)./evec5)
norm((hvec5-hvec6)./hvec5)


% construct reflected and transmitted field multipole expansions

[ra0E,rb0E]=rcoefs_die_arb12(nterms0,omega,radius0,eps0,cmu0,epsE,cmuE);
[ta0E,tb0E]=tcoefs_die_arb12(nterms0,omega,radius0,eps0,cmu0,epsE,cmuE);

[raa_diag,rbb_diag]=rmatr_diag(nterms0,ra0E,rb0E);
  atvec = em3sphlin(nterms0,aompole0);
  btvec = em3sphlin(nterms0,bompole0);
  aovec = atvec .* raa_diag;
  bovec = btvec .* rbb_diag;
  aimpoleE = em3linsph(nterms0,aovec);
  bimpoleE = em3linsph(nterms0,bovec);   

[taa_diag,tbb_diag]=rmatr_diag(nterms0,ta0E,tb0E);
  atvec = em3sphlin(nterms0,aompole0);
  btvec = em3sphlin(nterms0,bompole0);
  aovec = atvec .* taa_diag;
  bovec = btvec .* tbb_diag;
  aompoleE = em3linsph(nterms0,aovec);
  bompoleE = em3linsph(nterms0,bovec);   

'reflected field'
[evec7,hvec7]=em3taeval(rk0,center0,aimpoleE,bimpoleE,nterms0,target);
[evec7,hvec7]=em3dipole3ehimp(rk0,eps0,cmu0,evec7,hvec7)

'transmitted field'
[evec8,hvec8]=em3mpeval(rkE,center0,aompoleE,bompoleE,nterms0,target);
[evec8,hvec8]=em3dipole3ehimp(rkE,epsE,cmuE,evec8,hvec8)

% test jump condition, arbitrary target location on the sphere

rnorm = target-center;
rnorm = rnorm/norm(rnorm,2);

etangjump = cross(evec6+evec7,rnorm)-cross(evec8,rnorm)
htangjump = cross(hvec6+hvec7,rnorm)-cross(hvec8,rnorm)

enormjump = eps0*dot(evec6+evec7,rnorm)-epsE*dot(evec8,rnorm)
hnormjump = cmu0*dot(hvec6+hvec7,rnorm)-cmuE*dot(hvec8,rnorm)

