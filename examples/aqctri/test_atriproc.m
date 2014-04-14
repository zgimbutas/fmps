
%
%  Retrieve flat triangulation
%

geom_type = 2;
filename_geo = 'sphere80.a.tri';

[verts,ifaces,nverts,nfaces] = atriread(filename_geo);
nverts,nfaces

%
%  Construct triangle vertex, normal, area, and centroid arrays
%

ntri = nfaces;
[triangles,trianorm,triaarea,source]=atriproc(verts,ifaces);
