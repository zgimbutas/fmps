function plot_spheres(id,filename_geo,nspheres,center,radius)

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

trisurf(ifaces', verts(1,:)'*radius(i)+center(1,i), ...
       verts(2,:)'*radius(i)+center(2,i), verts(3,:)'*radius(i)+center(3,i),1)

end

axis equal
axis auto
alpha(.5)

figure(id)
hold off
