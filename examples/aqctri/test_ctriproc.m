
%
%  Retrieve cubic triangulation
%

geom_type = 4;
filename_geo = 'sphere80.c.tri';

[verts,ifaces,nverts,nfaces] = ctriread(filename_geo);
nverts,nfaces

%
%  Construct triangle vertex, normal, area, and centroid arrays
%

ntri = nfaces;
[triangles,trianorm,triaarea,source]=ctriproc(verts,ifaces);
