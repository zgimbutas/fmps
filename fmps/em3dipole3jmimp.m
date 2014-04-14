function [cjvec,cmvec]=em3dipole3jmimp(rk,eps,cmu,cjvec0,cmvec0)

cjvec=cjvec0*sqrt(cmu);
cjvec=cmvec0*sqrt(eps);
