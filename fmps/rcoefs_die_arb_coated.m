function [ra,rb,ta,tb]=rcoefs_die_arb_coated(nterms,omega,r1,r2,eps1,cmu1,eps2,cmu2,eps3,cmu3)
%RCOEFS_DIE_ARB_COATED Reflection coefficients for a coated sphere.
%
%  [ra,rb,ta,tb]=rcoefs_die_arb_coated ...
%           (nterms,omega,r1,r2,eps1,cmu1,eps2,eps2,eps3,cmu3);
%
%  Construct reflection and transmission coefficients
%  for an arbitrary coated dielectric sphere.
%
%  Incoming wave
%
%  Input parameters:
%
%    nterms - the number of terms in multipole expansions
%    omega - angular frequency parameter
%    r1 - inner sphere radius
%    r2 - outer sphere radius
%    eps1 - complex: permittivity of interior media (Layer 1)
%    cmu1 - complex: permeability of interior media (Layer 1)
%    eps2 - complex: permittivity of interior media (Layer 2)
%    cmu2 - complex: permeability of interior media (Layer 2)
%    eps3 - complex: permittivity of exterior media
%    cmu3 - complex: permeability of exterior media
%
%  Output parameters:
%
%    ra - complex(nterms+1): TM reflection coefficients
%    rb - complex(nterms+1): TE reflection coefficients
%    ta - complex(nterms+1): TM transmission coefficients
%    tb - complex(nterms+1): TE transmission coefficients
% 

[ra32,rb32] = rcoefs_die_arb21(nterms,omega,r2,eps3,cmu3,eps2,cmu2);
[ta32,tb32] = tcoefs_die_arb21(nterms,omega,r2,eps3,cmu3,eps2,cmu2);

[ra23,rb23] = rcoefs_die_arb12(nterms,omega,r2,eps2,cmu2,eps3,cmu3);
[ta23,tb23] = tcoefs_die_arb12(nterms,omega,r2,eps2,cmu2,eps3,cmu3);

[ra21,rb21] = rcoefs_die_arb21(nterms,omega,r1,eps2,cmu2,eps1,cmu1);

ra = ra32 + (ta23.*ra21.*ta32) ./ (1-ra23.*ra21);
rb = rb32 + (tb23.*rb21.*tb32) ./ (1-rb23.*rb21);

ta = ta21.*ta32 ./ (1-ra23.*ra21);
tb = tb21.*tb32 ./ (1-rb23.*rb21);
