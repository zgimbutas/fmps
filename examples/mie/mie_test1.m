% Test Mie scattering code
% Test extinction cross section calculation for PEC sphere.

omega=1;
eps0=1;
cmu0=1;
radius=2*pi;

nterms=100;
[a,b]=planew_ab(nterms);
[ra,rb]=rcoefs_pec(nterms,omega,radius,eps0,cmu0);

sa=a.*ra;
sb=b.*rb;

%%%q=sum( (2*(1:nterms)'+1).*real(sa+sb));
q=-real(a'*sa+b'*sb);

wavelength=2*pi/omega;
alpha=2*pi*radius/wavelength;

qext = q * (4/alpha^2);

