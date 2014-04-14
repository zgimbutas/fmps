
%
%  Retrieve cubic triangulation
%

geom_type = 4;
filename_geo = 'sphere80.c.tri';
%%filename_geo = 'sphere180.c.tri';
%%filename_geo = 'sphere720.c.tri';

[verts,ifaces,nverts,nfaces] = ctriread(filename_geo);
nverts,nfaces

%
%  Construct triangle vertex, normal, area, and centroid arrays
%

%ntri = nfaces;
%[triangles,trianorm,triaarea,source]=qtriproc(verts,ifaces);


%
%  Convert cubic triangulation in flat one, and store
%
[averts,iafaces]=ctri2atri(verts,ifaces,1);

geom_type = 2;
filename_geo_out = 'test1.a.tri';

atriwrite(filename_geo_out,averts,iafaces);
