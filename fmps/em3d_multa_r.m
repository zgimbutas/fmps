function [aompole,bompole,A]=em3d_multa_r(center,radius,atmpole,btmpole,nterms,A)

nspheres = A.nspheres;

ncoefs = (nterms+1)*(2*nterms+1);
aompole = zeros(ncoefs,nspheres);
bompole = zeros(ncoefs,nspheres);

for i=1:nspheres

if( A.type(i) == 1 ),

  [ra,rb] = rcoefs_pec(nterms,A.omega,radius(i),A.eps0,A.cmu0);
  [raa_diag,rbb_diag]=rmatr_diag(nterms,ra,rb);

  atvec = em3sphlin(nterms,atmpole(:,i));
  btvec = em3sphlin(nterms,btmpole(:,i));
  aovec = atvec .* raa_diag;
  bovec = btvec .* rbb_diag;
  aompole(:,i) = em3linsph(nterms,aovec);
  bompole(:,i) = em3linsph(nterms,bovec);   

end

if( A.type(i) == 2 ),

  [ra,rb] = rcoefs_die_arb(nterms,A.omega,radius(i),A.eps0,A.cmu0,...
            A.eps(i),A.cmu(i));
  [raa_diag,rbb_diag]=rmatr_diag(nterms,ra,rb);

  atvec = em3sphlin(nterms,atmpole(:,i));
  btvec = em3sphlin(nterms,btmpole(:,i));
  aovec = atvec .* raa_diag;
  bovec = btvec .* rbb_diag;
  aompole(:,i) = em3linsph(nterms,aovec);
  bompole(:,i) = em3linsph(nterms,bovec);

end

if( A.type(i) == 3 ),

  rr = (nterms+1)^2;

  ampole = emabrotaf(nterms,A.rot(:,i),atmpole(:,i));
  bmpole = emabrotaf(nterms,A.rot(:,i),btmpole(:,i));
  atvec = em3sphlin(nterms,ampole);
  btvec = em3sphlin(nterms,bmpole);
  aovec = A.raa(1:rr,1:rr) * atvec + A.rab(1:rr,1:rr) * btvec;
  bovec = A.rba(1:rr,1:rr) * atvec + A.rbb(1:rr,1:rr) * btvec;
  ampole = em3linsph(nterms,aovec);
  bmpole = em3linsph(nterms,bovec);
  aompole(:,i) = emabrotab(nterms,A.rot(:,i),ampole);
  bompole(:,i) = emabrotab(nterms,A.rot(:,i),bmpole);
  
end

end

