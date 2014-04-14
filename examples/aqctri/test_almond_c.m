
%
%  Retrieve cubic triangulation
%

geom_type = 4;
filename_geo = 'sphere320.c.tri';
%filename_geo = 'sphere1280.c.tri';
%filename_geo = 'sphere5120.c.tri';
%filename_geo = 'sphere20480.c.tri';

[verts,ifaces,nverts,nfaces] = ctriread(filename_geo);
nverts,nfaces


%
%  Construct NASA almond
%
verts=sphere2almond(verts);

%
%  Store cubic triangulation
%
geom_type = 4;
filename_geo_out = 'almond0.c.tri';

ctriwrite(filename_geo_out,verts,ifaces);


%
%  Convert cubic triangulation in flat one, and store
%
[averts,iafaces]=ctri2atri(verts,ifaces,1);

geom_type = 2;
filename_geo_out = 'test0.a.tri';

atriwrite(filename_geo_out,averts,iafaces);
