
for itype=1:4,

[verts,ifaces,nverts,nfaces] = rsolid(itype);
nverts,nfaces

geom_type = 2;
filename_geo_out = ['test', num2str(itype), '.a.tri'];

atriwrite(filename_geo_out,verts,ifaces);

end
