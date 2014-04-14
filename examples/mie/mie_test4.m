%
%  Mie scattering from an arbitrary dielectric sphere
%  gold at 620nm, 500nm, plasmon resonance, iron

omega=1;
eps0=1;
cmu0=1;

%%eps1=(0.2039+3.3056i)^2; % gold at 620nm 
eps1=(0.9711+1.8563i)^2; % gold at 500nm 
%%eps1=-2+.01i; % plasmon resonance
%%eps1=(1.27+1.37i)^2; % iron
cmu1=1;

n=2000;
alphas=zeros(n,1);
qext=zeros(n,1);
qsca=zeros(n,1);

nterms=100;
[a,b]=planew_ab(nterms);

for i=1:n

radius=(i/n) * 1*pi;
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


qabs=qext-qsca;

plot(alphas/pi,qsca,'+',alphas/pi,qext,'*',alphas/pi,qabs,'.');
