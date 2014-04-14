%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  STEP 3. Read in wavelength data from gold_jc_data and initiallize
%          output arrays for scattering absorption, extinction as well as
%          electric and magnetic dipole moments.
%         

gold_jc_data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  STEP 4. Run solver (using gold_jc_data points) for a single frequency
%

ifreq = 50;

wavelength=gold_jc(ifreq,2)
omega=2*pi/wavelength

epsE=1;
cmuE=1;
rkE=omega*sqrt(epsE)*sqrt(cmuE);

eps0=1.2;
cmu0=1;
rk0=omega*sqrt(eps0)*sqrt(cmu0);

center0=[0,0,0]';
radius0=200;

ntermsE=ceil(abs(radius0*rkE)*1.2)+24*3
nterms0=ceil(abs(radius0*rk0)*1.2)+24*3


re_n=gold_jc(ifreq,3);
im_n=gold_jc(ifreq,4);

eps1=(re_n+ima*im_n)^2
cmu1=1

r1=50;
r2=200;


% construct reflected and transmitted field multipole expansions

% reflected
[aimpoleE,bimpoleE]=em3formta(rkE,source,cjvec,cmvec,center0,nterms0);

[raE0,rbE0]=rcoefs_die_arb21(ntermsE,omega,radius0,epsE,cmuE,eps0,cmu0);

[raa_diag,rbb_diag]=rmatr_diag(ntermsE,raE0,rbE0);
  atvec = em3sphlin(ntermsE,aimpoleE);
  btvec = em3sphlin(ntermsE,bimpoleE);
  aovec = atvec .* raa_diag;
  bovec = btvec .* rbb_diag;
  aompoleE = em3linsph(ntermsE,aovec);
  bompoleE = em3linsph(ntermsE,bovec);   



[ra,rb,ta,tb]=rcoefs_die_arb_coated(nterms,omega,r1,r2,eps1,cmu1,eps0,cmu0,epsE,cmuE);



ai=reshape(aimpoleE,nterms0+1,2*nterms0+1);
ao=reshape(aompoleE+aompole0E,nterms0+1,2*nterms0+1);

'ra, approx'
%ao(2,nterms0:nterms0+2) ./  ai(2,nterms0:nterms0+2)
%ao(3,nterms0:nterms0+2) ./  ai(3,nterms0:nterms0+2)
%ao(4,nterms0:nterms0+2) ./  ai(4,nterms0:nterms0+2)
ao(2,nterms0+2) ./  ai(2,nterms0+2)
ao(3,nterms0+2) ./  ai(3,nterms0+2)
ao(4,nterms0+2) ./  ai(4,nterms0+2)

'ra, exact'
ra(1:3)


if( 2 == 2 ), 
bi=reshape(bimpoleE,nterms0+1,2*nterms0+1);
bo=reshape(bompoleE+bompole0E,nterms0+1,2*nterms0+1);

'rb, approx'
%bo(2,nterms0:nterms0+2) ./  bi(2,nterms0:nterms0+2)
%bo(3,nterms0:nterms0+2) ./  bi(3,nterms0:nterms0+2)
%bo(4,nterms0:nterms0+2) ./  bi(4,nterms0:nterms0+2)
bo(2,nterms0+2) ./  bi(2,nterms0+2)
bo(3,nterms0+2) ./  bi(3,nterms0+2)
bo(4,nterms0+2) ./  bi(4,nterms0+2)
end

'rb, exact'
rb(1:3)

