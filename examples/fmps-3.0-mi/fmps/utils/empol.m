function [pvec,mvec]=empol(wavelength,rk,ampole,bmpole,nterms,radius)
%EMPOL: Evaluate electric polarization and magnetization vectors.
%
%  Evaluate electric polarization and magnetization vectors 
%  for a single sphere.
%
%  [pvec,mvec]=empol(wavelength,rk,ampole,bmpole,nterms,radius);
%
%  Input parameters:
%
%    wavelength - unused
%    rk - the frequency parameter
%    ampole,bmpole - complex(ncoefs): outgoing EM multipoles
%    nterms - the number of terms in multipole expansion
%    radius - unused
%
%  Output parameter:
%
%    pvec - complex(3) - electric polarization vector
%    mvec - complex(3) - magnetic polarization (magnetization) vector
%
%  All EM multipoles are NCOEFS = (NTERMS+1)-by-(2*NTERMS+1) complex matrices.
%

pvec = zeros(3,1) + 1i*zeros(3,1);
mvec = zeros(3,1) + 1i*zeros(3,1);

mex_id_ = 'empol(i double[x], i dcomplex[x], i dcomplex[], i dcomplex[], i int[x], i double[x], io dcomplex[], io dcomplex[])';
[pvec, mvec] = fmps_r2012a(mex_id_, wavelength, rk, ampole, bmpole, nterms, radius, pvec, mvec, 1, 1, 1, 1);



