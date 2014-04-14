%
%  SET ALL PARAMETERS
%

B.itype=1;
B.nquad=min(12,nterms);
B.nphi=2*B.nquad+1;
B.ntheta=B.nquad+1;

[B.rnodes,B.weights,B.nnodes]=e3fgrid(B.itype,B.nquad,B.nphi,B.ntheta);

x = reshape(B.rnodes(1,:),B.nphi,B.ntheta);
y = reshape(B.rnodes(2,:),B.nphi,B.ntheta);
z = reshape(B.rnodes(3,:),B.nphi,B.ntheta);

%%%surf(x,y,z)
%figure(1)
%surf([x; x(1,:)],[y; y(1,:)],[z; z(1,:)]);
%axis equal

xp = [x; x(1,:)];
yp = [y; y(1,:)];
zp = [z; z(1,:)];


%
%  Evaluate E and H fields on all spheres
%
eivecs = zeros(3,B.nnodes,nspheres) + 1i*zeros(3,B.nnodes,nspheres);
hivecs = zeros(3,B.nnodes,nspheres) + 1i*zeros(3,B.nnodes,nspheres);
esvecs = zeros(3,B.nnodes,nspheres) + 1i*zeros(3,B.nnodes,nspheres);
hsvecs = zeros(3,B.nnodes,nspheres) + 1i*zeros(3,B.nnodes,nspheres);
eovecs = zeros(3,B.nnodes,nspheres) + 1i*zeros(3,B.nnodes,nspheres);
hovecs = zeros(3,B.nnodes,nspheres) + 1i*zeros(3,B.nnodes,nspheres);

for i=1:nspheres

 [evecs,hvecs]=...
    em3taevaleh(zk,center(:,i),aimpole(:,i),bimpole(:,i),nterms, ...
    center(:,i),radius(i),B.rnodes,B.weights,B.nphi,B.ntheta);
  eivecs(:,:,i) = evecs;
  hivecs(:,:,i) = hvecs;

 [evecs,hvecs]=...
    em3taevaleh(zk,center(:,i),asmpole(:,i),bsmpole(:,i),nterms, ...
    center(:,i),radius(i),B.rnodes,B.weights,B.nphi,B.ntheta);
  esvecs(:,:,i) = evecs;
  hsvecs(:,:,i) = hvecs;

 [evecs,hvecs]=...
    em3mpevaleh(zk,center(:,i),aompole(:,i),bompole(:,i),nterms, ...
    center(:,i),radius(i),B.rnodes,B.weights,B.nphi,B.ntheta);
  eovecs(:,:,i) = evecs;
  hovecs(:,:,i) = hvecs;

end


eivecs = reshape(eivecs,3,B.nphi,B.ntheta,nspheres);
hivecs = reshape(hivecs,3,B.nphi,B.ntheta,nspheres);
esvecs = reshape(esvecs,3,B.nphi,B.ntheta,nspheres);
hsvecs = reshape(hsvecs,3,B.nphi,B.ntheta,nspheres);
eovecs = reshape(eovecs,3,B.nphi,B.ntheta,nspheres);
hovecs = reshape(hovecs,3,B.nphi,B.ntheta,nspheres);

% total field 
etvecs = eivecs + esvecs + eovecs;
htvecs = hivecs + hsvecs + hovecs;

% scattered field
eavecs = esvecs + eovecs;
havecs = hsvecs + hovecs;

xp = [x; x(1,:)];
yp = [y; y(1,:)];
zp = [z; z(1,:)];


%
%  Plot Re(E_y^{incoming}) on all spheres
%
figure(101);
for isphere = 1:nspheres
    
f = squeeze((eivecs(2,:,:,isphere)));
fp = [f; f(1,:)];

xc = xp*radius(isphere) + center(1,isphere);
yc = yp*radius(isphere) + center(2,isphere);
zc = zp*radius(isphere) + center(3,isphere);

surf(xc,yc,zc,real(fp));
axis equal
hold on

end
colorbar
xlabel('x')
ylabel('y')
zlabel('z')
title('Re(E_y^{incoming})')
figure(101)


%
%  Plot Re(E_y^{total}) on all spheres
%
figure(102);
for isphere = 1:nspheres
    
f = squeeze((etvecs(2,:,:,isphere)));
fp = [f; f(1,:)];

xc = xp*radius(isphere) + center(1,isphere);
yc = yp*radius(isphere) + center(2,isphere);
zc = zp*radius(isphere) + center(3,isphere);

surf(xc,yc,zc,real(fp));
axis equal
hold on

end
colorbar
xlabel('x')
ylabel('y')
zlabel('z')
title('Re(E_y^{total})')
figure(102)


%
%  Plot Re(E_y^{scat}) on all spheres
%
figure(103);
for isphere = 1:nspheres
    
f = squeeze((esvecs(2,:,:,isphere)));
fp = [f; f(1,:)];

xc = xp*radius(isphere) + center(1,isphere);
yc = yp*radius(isphere) + center(2,isphere);
zc = zp*radius(isphere) + center(3,isphere);

surf(xc,yc,zc,real(fp));
axis equal
hold on

end
colorbar
xlabel('x')
ylabel('y')
zlabel('z')
title('Re(E_y^{scat})')
figure(103)

