function [rnodes,weights,nnodes]=e3fgrid(itype,nquad,nphi,ntheta)
%E3FGRID: Construct a spherical grid in R^3.
%
%  [rnodes,weights,nnodes]=e3fgrid(itype,nquad,nphi,ntheta);
%
%  Construct a spherical grid on a unit sphere in R^3. For itype=1, the nodes
%  in theta direction are Gauss-Legendre quadrature nodes in z, and uniformly
%  spaced in phi (compatible with FFTPACK).
%
%  Input parameters:
%
%    itype - (internal parameter, leave set to 1).
%    nquad - (internal parameter, leave set to ntheta).
%    nphi - number of nodes in angular phi direction
%    ntheta - number of nodes in theta direction 
%
%  Output parameters:
%
%    rnodes - real (3,nphi*ntheta): the spherical grid
%    weights - real (nphi*ntheta): the associated quadrature weights
%    nnodes - the number of nodes in spherical grid
%    

rnodes = zeros(3,nphi*ntheta);
weights = zeros(1,nphi*ntheta);
nnodes = 0;

mex_id_ = 'e3fgrid(i int[x], i int[x], i int[x], i int[x], i double[], io double[], io int[])';
[weights, nnodes] = fmps(mex_id_, itype, nquad, nphi, ntheta, rnodes, weights, nnodes, 1, 1, 1, 1);


