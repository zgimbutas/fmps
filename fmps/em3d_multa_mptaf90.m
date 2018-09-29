function [asmpole,bsmpole] = em3d_multa_mptaf90(nspheres,nterms,ncoefs,omega,eps0,cmu0,center,radius,aompole,bompole,rnodes,weights,nphi,ntheta)
%EM3D_MULTA_MPTAF90: Direct (slow) EM multipole to local translation routine.
%
%  [ASMPOLE,BSMPOLE] = EM3D_MULTA_MPTAF90(NSPHERES,NTERMS,NCOEFS,...
%       OMEGA,EPS0,CMU0,CENTER,RADIUS,AOMPOLE,BOMPOLE,...
%       RNODES,WEIGHTS,NPHI,NTHETA);
%
%  Convert outgoing EM multipole expansions for a collection of spheres to
%  the incoming EM multipole expansions. Self interactions are NOT included.
%  
%  All EM multipoles are (NTERMS+1)-by-(2*NTERMS+1) complex matrices.
%
%  Input parameters:
%
%    nspheres - the number of spheres
%    nterms - the number of terms in multipole expansions
%    ncoefs - currently, must be set to to (NTERMS+1)*(2*NTERMS+1)
%    omega - angular frequency parameter
%    eps0 - complex: permittivity of exterior media
%    cmu0 - complex: permeability of exterior media
%    center - real(3,nsource): sphere center locations
%    radius - real(nsource): sphere radii
%    aompole,bompole - complex(ncoefs,nspheres): outgoing EM multipoles
%    rnodes,weights,nphi,ntheta - 
%          spherical grid, constructed via a preceding call to e3fgrid(nterms)
%                 
%  Output parameters:
%
%    asmpole,bsmpole - complex(ncoefs,nspheres): incoming EM multipoles
%
%

asmpole = zeros(ncoefs,nspheres) + 1i*zeros(ncoefs,nspheres);
bsmpole = zeros(ncoefs,nspheres) + 1i*zeros(ncoefs,nspheres);

if( nspheres > 1 ),
mex_id_ = 'em3dmpta(i int[x], i int[x], i int[x], i double[x], i dcomplex[x], i dcomplex[x], i double[], i double[], i dcomplex[], i dcomplex[], io dcomplex[], io dcomplex[], i double[], i double[], i int[x], i int[x])';
[asmpole, bsmpole] = fmps(mex_id_, nspheres, nterms, ncoefs, omega, eps0, cmu0, center, radius, aompole, bompole, asmpole, bsmpole, rnodes, weights, nphi, ntheta, 1, 1, 1, 1, 1, 1, 1, 1);
end


