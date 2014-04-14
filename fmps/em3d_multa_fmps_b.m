function [y,A]=em3d_multa_fmps_b(A,x)
%
%  EM multisphere solver multiplication routine
%
%  (E,H)^{scat}_{self} - R MPTA (E,H)^{scat}_{self} - R TATA R_big MPMP (E,H)^{scat}_{self}
%
'em3d_multa_fmps_b'

x0 = reshape(x,A.ncoefs,A.nspheres,2);

aimpole = x0(:,:,1);
bimpole = x0(:,:,2);

imulta_type = 3;
t1=toc;
if( imulta_type == 1 ),
% Matlab call, very slow...
[asmpole,bsmpole]=em3d_multa_mpta(A.center,A.radius,aimpole,bimpole,A.nterms,A);
end

if( imulta_type == 2 ),
% Fortran 90, same calculation
[asmpole,bsmpole] = em3d_multa_mptaf90(A.nspheres,A.nterms,A.ncoefs,A.omega,A.eps0,A.cmu0,A.center,A.radius,aimpole,bimpole,A.rnodes,A.weights,A.nphi,A.ntheta);
end

if( imulta_type == 3 ),
% Fortran 90, same calculation with FMM
[asmpole,bsmpole] = em3d_multa_mptafmm(A.iprec,A.nspheres,A.nterms,A.ncoefs,A.omega,A.eps0,A.cmu0,A.center,A.radius,aimpole,bimpole,A.rnodes,A.weights,A.nphi,A.ntheta);
end
'fmm',toc-t1

% Big sphere contribution

omega=A.omega;
rkE=A.rkE;
rk0=A.rk0;
epsE=A.epsE;
cmuE=A.cmuE;
eps0=A.eps0;
cmu0=A.cmu0;
ncoefs0=A.ncoefs0;
nterms0=A.nterms0;
center0=A.center0;
radius0=A.radius0;
nspheres=A.nspheres;
nterms=A.nterms;

% shift outgoing expansions to the exterior sphere
aompole0 = zeros(ncoefs0,1);
bompole0 = zeros(ncoefs0,1);

% truncate outgoing multipole expansions on the big sphere, 
% this is much faster if than naive translation if nterms0 is big
  t1=toc;
if_truncate = 1;

if( if_truncate == 0 ),

G0=get_e3fgrid(nterms0);

for i=1:nspheres

[ampoletmp,bmpoletmp]=...
      em3mpmp3(rk0,A.center(:,i),aimpole(:,i),bimpole(:,i),A.nterms, ...
      center0,nterms0,radius0,...
      G0.rnodes,G0.weights,G0.nphi,G0.ntheta);

aompole0 = aompole0+ampoletmp;
bompole0 = bompole0+bmpoletmp;

end

end

if( if_truncate == 1 ),

G0=get_e3fgrid_trunc(nterms0,nterms);

for i=1:nspheres

[ampoletmp,bmpoletmp]=...
      em3mpmp3_trunc(rk0,A.center(:,i),aimpole(:,i),bimpole(:,i),A.nterms, ...
      center0,nterms0,radius0,...
      G0.rnodes,G0.weights,G0.nphi,G0.ntheta);

aompole0 = aompole0+ampoletmp;
bompole0 = bompole0+bmpoletmp;

end

end

'mpmp',toc-t1
  t1=toc;
% construct the reflected field multipole expansion
[ra0E,rb0E]=rcoefs_die_arb12(nterms0,omega,radius0,eps0,cmu0,epsE,cmuE);

[raa_diag,rbb_diag]=rmatr_diag(nterms0,ra0E,rb0E);
  atvec = em3sphlin(nterms0,aompole0);
  btvec = em3sphlin(nterms0,bompole0);
  aovec = atvec .* raa_diag;
  bovec = btvec .* rbb_diag;
  asmpole0E = em3linsph(nterms0,aovec);
  bsmpole0E = em3linsph(nterms0,bovec);   

'r_big',toc-t1
% shift the incoming expansion to interior spheres
  t1=toc;
if( if_truncate == 0 ),
G=get_e3fgrid(nterms);

for i=1:nspheres

[ampoletmp,bmpoletmp]=...
      em3tata3(rk0,center0,asmpole0E,bsmpole0E,nterms0, ...
      A.center(:,i),A.nterms,A.radius(i),...
      G.rnodes,G.weights,G.nphi,G.ntheta);

asmpole(:,i) = asmpole(:,i)+ampoletmp;
bsmpole(:,i) = bsmpole(:,i)+bmpoletmp;

end
end

if( if_truncate == 1 ),
G=get_e3fgrid(nterms);

for i=1:nspheres

[ampoletmp,bmpoletmp]=...
      em3tata3_trunc(rk0,center0,asmpole0E,bsmpole0E,nterms0, ...
      A.center(:,i),A.nterms,A.radius(i),...
      G.rnodes,G.weights,G.nphi,G.ntheta);

asmpole(:,i) = asmpole(:,i)+ampoletmp;
bsmpole(:,i) = bsmpole(:,i)+bmpoletmp;

end
end

'tata',toc-t1

  t1=toc;
% Apply reflection matrices
[aompole,bompole,A]=em3d_multa_r(A.center,A.radius,asmpole,bsmpole,A.nterms,A);
'r_small',toc-t1
atmpole = aimpole - aompole;
btmpole = bimpole - bompole;

y = reshape([atmpole btmpole], A.ncoefs*A.nspheres*2, 1);

toc
