
itype=4;

[verts,ifaces,nverts,nfaces] = rsolid(itype);
nverts,nfaces

geom_type = 2;
filename_geo_out = ['test_sphere', num2str(nfaces), '.a.tri'];

atriwrite(filename_geo_out,verts,ifaces);


for i=1:4

[verts1,ifaces1]=atrirefine(verts,ifaces);
verts = verts1;
ifaces = ifaces1;

% map to the unit sphere
verts = verts ./ sqrt(repmat(sum(verts.^2, 1),3,1));

nverts = size(verts,2)
nfaces = size(ifaces,2)

geom_type = 2;
filename_geo_out = ['test_sphere', num2str(nfaces), '.a.tri'];

atriwrite(filename_geo_out,verts,ifaces);

end