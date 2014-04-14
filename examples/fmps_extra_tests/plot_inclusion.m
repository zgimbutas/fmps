%
%  SET ALL PARAMETERS
%

B.itype=1;
B.nquad=12;
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
%  Plot all spheres
%
filename_geo = '2ellipsoids-25x25x75-sep40-draft.a.tri'
[verts,ifaces,nverts,nfaces] = atriread(filename_geo);

figure(4);
for isphere = 1:nspheres
    
xc = xp*radius(isphere) + center(1,isphere);
yc = yp*radius(isphere) + center(2,isphere);
zc = zp*radius(isphere) + center(3,isphere);

figure(4)
%surf(xc,yc,zc,ones(B.nphi+1,B.ntheta));
hold on
axis equal

gamma0 = A.rot(3,isphere);
beta0 = A.rot(2,isphere);
alpha0 = A.rot(1,isphere);

rot1 = [cos(gamma0) sin(gamma0) 0;
       -sin(gamma0) cos(gamma0) 0;
        0          0          1];
rot2 = [cos(beta0) 0 sin(beta0);
        0         1        0 ;
       -sin(beta0) 0 cos(beta0)];
rot3 = [cos(alpha0) sin(alpha0) 0;
       -sin(alpha0) cos(alpha0) 0;
        0          0          1];
rot = rot3 * rot2 * rot1;

verts_rot = rot * verts;

figure(4)
trisurf(ifaces', (center(1,isphere)+verts_rot(1,:))', ...
      (center(2,isphere)+verts_rot(2,:))', ...
      (center(3,isphere)+verts_rot(3,:))', ...
      ones(nfaces,1))
hold on
axis equal

end
title('Geometry')
xlabel('x')
ylabel('y')
zlabel('z')
%alpha(0.5)
figure(4)

