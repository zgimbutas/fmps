function [evec,hvec]=em3dipole3mtimp(rk,eps,cmu,source,target,cmvec)

cmvec=cmvec*sqrt(eps);
[evec,hvec] = em3dipole3mt(rk,source,target,cmvec);
evec=evec/sqrt(eps);
hvec=hvec/sqrt(cmu);
