%
%  Mie scattering from an arbitrary dielectric sphere
%  Coated gold at 500nm
%  NaNs for small spheres (needs to be fixed)

omega=1;
eps0=1;
cmu0=1;

eps1=1.33^2; % water
cmu1=1;

eps2=(0.9711+1.8563i)^2; % gold at 500nm 
cmu2=1;

n=2000;
alphas=zeros(n,1);
qext=zeros(n,1);
qsca=zeros(n,1);

nterms=100;
[a,b]=planew_ab(nterms);

for i=1:n

radius1=(i/n) * 1*pi *0.80 *10;
radius2=(i/n) * 1*pi *1.00 *10;
wavelength=2*pi/omega;
alpha=2*pi*radius2/wavelength;

[ra,rb]=rcoefs_die_arb_coated(nterms,omega,radius1,radius2,eps2,cmu2,eps1,cmu1,eps0,cmu0);

k=find(isnan(ra)); ra(k) = 0;
k=find(isnan(rb)); rb(k) = 0;

sa=a.*ra;
sb=b.*rb;

q=-real(a'*sa+b'*sb);

qext(i) = q *(4/alpha^2);

q=sum((abs(sa).^2+abs(sb).^2));

qsca(i) = q *(4/alpha^2);

alphas(i) = alpha;

end


qabs=qext-qsca;

plot(alphas/pi,qsca,'+',alphas/pi,qext,'*',alphas/pi,qabs,'.');
