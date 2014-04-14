
%
%  Retrieve quadratic triangulation
%

geom_type = 3;
filename_geo = 'sphere80.q.tri';

[verts,ifaces,nverts,nfaces] = qtriread(filename_geo);
nverts,nfaces

%
%  Construct triangle vertex, normal, area, and centroid arrays
%

ntri = nfaces;
[triangles,trianorm,triaarea,source]=qtriproc(verts,ifaces);
