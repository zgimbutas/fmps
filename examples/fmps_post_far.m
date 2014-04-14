%
%  SET ALL PARAMETERS
%

radius_far = max(sqrt(sum(center.^2,1)) + radius)
rk0=omega*sqrt(eps0)*sqrt(cmu0)
nquad_far=ceil(abs(radius_far*rk0)*1.2)+36


C.itype=1;
C.nquad=nquad_far;
C.nphi=2*C.nquad+1;
C.ntheta=C.nquad+1;

[C.rnodes,C.weights,C.nnodes]=e3fgrid(C.itype,C.nquad,C.nphi,C.ntheta);


efar = zeros(3,C.nnodes) + 1i*zeros(3,C.nnodes);
hfar = zeros(3,C.nnodes) + 1i*zeros(3,C.nnodes);

for i=1:nspheres

  [evecs,hvecs]=...
    em3mpfareh(rk0,center(:,i),aompole(:,i),bompole(:,i),nterms, ...
    C.rnodes,C.weights,C.nphi,C.ntheta);

  efar = efar + evecs;
  hfar = hfar + hvecs;

end




%
%  Plot abs(E_{far}) 
%
figure(4);
    
f = efar;
f = sqrt(sum(abs(f).^2,1));

x = reshape(C.rnodes(1,:),C.nphi,C.ntheta);
y = reshape(C.rnodes(2,:),C.nphi,C.ntheta);
z = reshape(C.rnodes(3,:),C.nphi,C.ntheta);
f = reshape(f,C.nphi,C.ntheta);

xp = [x; x(1,:)] *radius_far;
yp = [y; y(1,:)] *radius_far;
zp = [z; z(1,:)] *radius_far;
fp = [f; f(1,:)];

surf(xp,yp,zp,abs(fp));
axis equal

colorbar
xlabel('x')
ylabel('y')
zlabel('z')
title('abs(E_{far})')
figure(4)




%
%  Plot 20*log10(abs(E_{far}))
%
figure(6);

surf(xp,yp,zp,20*log10(abs(fp)));
axis equal

colorbar
xlabel('x')
ylabel('y')
zlabel('z')
title('20*log10(abs(E_{far}))')
figure(6)



