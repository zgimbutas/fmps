%
% Test EM multipole to local translation operator
%

cjvec=[1;1;1];
cmvec=[0;0;0];

center=[0;0;0];
source=center+[0.1;0.2;-0.3];

center1=[10;-20;30];
target=center1+[.1;-.2;.3];

rk=1;
nterms=12;

[ampole,bmpole]=em3formmp(rk,source,cjvec,cmvec,center,nterms);

'directly'
[evec1,hvec1]=dipole3et(rk,source,target,cjvec)

'via mpole'
[evec2,hvec2]=em3mpeval(rk,center,ampole,bmpole,nterms,target)

'errors'
evec1-evec2
hvec1-hvec2


%
% Apply translation operator
%

itype=1;
nquad=nterms;
nphi=2*nquad+1;
ntheta=nquad+1;

A.nterms = nterms;
A.nquad = nquad;
A.nphi = nphi;
A.ntheta = ntheta;

[A.rnodes,A.weights,A.nnodes]=e3fgrid(itype,nquad,nphi,ntheta);

radius1=3; nterms1=12;
[ampole1,bmpole1]=em3mpta3(rk,center,ampole,bmpole,nterms,center1,nterms1,radius1,A.rnodes,A.weights,nphi,ntheta);



'via local'
[evec3,hvec3]=em3taeval(rk,center1,ampole1,bmpole1,nterms1,target)

'errors'
evec1-evec3
hvec1-hvec3

