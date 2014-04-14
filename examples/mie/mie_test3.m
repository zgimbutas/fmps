%
%  Mie scattering from an arbitrary dielectric sphere
%  Non-absorbing modest constrast ratio (water, n=1.33)

omega=1;
eps0=1;
cmu0=1;

eps1=1.33^2; % water
cmu1=1;

n=1000;
alphas=zeros(n,1);
qext=zeros(n,1);
qsca=zeros(n,1);

nterms=100;
[a,b]=planew_ab(nterms);

for i=1:n

radius=(i/n) * 10*pi;
wavelength=2*pi/omega;
alpha=2*pi*radius/wavelength;

[ra,rb]=rcoefs_die_arb(nterms,omega,radius,eps0,cmu0,eps1,cmu1);

sa=a.*ra;
sb=b.*rb;

q=-real(a'*sa+b'*sb);

qext(i) = q *(4/alpha^2);

q=sum((abs(sa).^2+abs(sb).^2));

qsca(i) = q *(4/alpha^2);

alphas(i) = alpha;

end


plot(alphas/pi,qsca,'+');
%%plot(alphas/pi,qext,'+');
