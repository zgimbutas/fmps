function [evec,hvec]=dipole3mtimp(rk,eps,cmu,source,target,cmvec)

cmvec=cmvec*sqrt(eps);
[evec,hvec] = dipole3mt(rk,source,target,cmvec);
evec=evec/sqrt(eps);
hvec=hvec/sqrt(cmu);
