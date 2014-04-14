function A=get_e3fgrid(nterms)
%GET_E3FGRID: Construct a spherical grid structure in R^3.
%
%  A=get_e3fgrid(nterms);
%
%  Construct a spherical grid on a unit sphere in R^3. For itype=1, the nodes
%  in theta direction are Gauss-Legendre quadrature nodes in z, and uniformly
%  spaced in phi (compatible with FFTPACK).
%
%  Input parameters:
%
%    nterms - number of terms
%
%  Output parameters:
%
%    A - data structure
%    A.nphi - number of terms in angular phi direction (2*nterms+1)
%    A.ntheta - number of nodes in theta direction (nterms+1)
%    A.rnodes - real (3,nphi*ntheta): the spherical grid
%    A.weights - real (nphi*ntheta): the associated quadrature weights
%    A.nnodes - the number of nodes in spherical grid
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

