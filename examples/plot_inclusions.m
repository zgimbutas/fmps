function plot_inclusions(id,filename_geo,nspheres,center,radius,sphere_rot)

%
%  Retrieve flat triangulation
%

geom_type = 2;

fid = fopen(filename_geo,'r');

nverts=0;
nfaces=0;
[nverts] = fscanf(fid,'%d',1);
[nfaces] = fscanf(fid,'%d',1);

verts=zeros(3,nverts);
ifaces=zeros(3,nfaces);

[verts] = fscanf(fid,'%f',[3,nverts]);
[ifaces] = fscanf(fid,'%d',[3,nfaces]);

fclose(fid);


figure(id)
hold on

for i=1:nspheres

rotmat = em3orient(sphere_rot(:,i),1);
verts_rot = rotmat*verts;

trisurf(ifaces', verts_rot(1,:)'+center(1,i), ...
        verts_rot(2,:)'+center(2,i), verts_rot(3,:)'+center(3,i),1)

end

axis equal
axis auto
hold off
xlabel('x')
ylabel('y')
zlabel('z')
figure(id)
