function [evec,hvec]=em3dipole3etimp(rk,eps,cmu,source,target,cjvec)

cjvec=cjvec*sqrt(cmu);
[evec,hvec] = em3dipole3et(rk,source,target,cjvec);
evec=evec/sqrt(eps);
hvec=hvec/sqrt(cmu);
