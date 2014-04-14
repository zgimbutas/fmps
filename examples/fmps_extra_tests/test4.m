%
% Test EM local expansion
%

cjvec=[1;1;1];
cmvec=[0;0;0];

radius=300;
target=[0;-radius;0];
center=[0;0;0];

source=[1000;-2000;3000];

omega=1;
eps=1;
cmu=1

rk=omega*sqrt(eps)*sqrt(cmu);

%%%nterms=1120;
nterms=ceil(abs(radius*rk)*1.2);


[ampole,bmpole]=em3formta(rk,source,cjvec,cmvec,center,nterms);

'directly'
[evec1,hvec1]=dipole3et(rk,source,target,cjvec)

'via local'
[evec2,hvec2]=em3taeval(rk,center,ampole,bmpole,nterms,target)

'errors'
evec1-evec2
hvec1-hvec2

'relative errors'
norm((evec1-evec2)./evec1)
norm((hvec1-hvec2)./hvec1)


% Reshape multipole arrays as (NTERMS+1)-by-(2*NTERMS+1) complex matrices.
amatr=reshape(ampole,nterms+1,2*nterms+1);
bmatr=reshape(bmpole,nterms+1,2*nterms+1);
