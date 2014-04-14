function [evec,hvec]=em3dipole3ehimp(rk,eps,cmu,evec0,hvec0)

evec=evec0/sqrt(eps);
hvec=hvec0/sqrt(cmu);

