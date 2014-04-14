%
%  Mie scattering from PEC sphere
%

omega=1;
eps0=1;
cmu0=1;

n=300;
alphas=zeros(n,1);
qext=zeros(n,1);
qsca=zeros(n,1);

nterms=100;
[a,b]=planew_ab(nterms);

for i=1:n

radius=i/n * 4*pi;
wavelength=2*pi/omega;
alpha=2*pi*radius/wavelength;

[ra,rb]=rcoefs_pec(nterms,omega,radius,eps0,cmu0);

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
