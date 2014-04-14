%
% Test EM local expansion
%

cjvec=[1;1;1];
cmvec=[0;0;0];

target=[0.1;0.2;-0.3];
center=[0;0;0];

source=[10;-20;30];

rk=1;
nterms=12;

[ampole,bmpole]=em3formta(rk,source,cjvec,cmvec,center,nterms);

'directly'
[evec1,hvec1]=dipole3et(rk,source,target,cjvec)

'via local'
[evec2,hvec2]=em3taeval(rk,center,ampole,bmpole,nterms,target)

'errors'
evec1-evec2
hvec1-hvec2


% Reshape multipole arrays as (NTERMS+1)-by-(2*NTERMS+1) complex matrices.
amatr=reshape(ampole,nterms+1,2*nterms+1);
bmatr=reshape(bmpole,nterms+1,2*nterms+1);
