
%
%  Retrieve quadratic triangulation
%

geom_type = 3;
filename_geo = 'sphere320.q.tri';
%filename_geo = 'sphere1280.q.tri';
%filename_geo = 'sphere5120.q.tri';
%filename_geo = 'sphere20480.q.tri';

[verts,ifaces,nverts,nfaces] = qtriread(filename_geo);
nverts,nfaces


%
%  Construct NASA almond
%
verts=sphere2almond(verts);

%
%  Store quadratic triangulation
%
geom_type = 3;
filename_geo_out = 'almond0.q.tri';

qtriwrite(filename_geo_out,verts,ifaces);


%
%  Convert quadratic triangulation in flat one, and store
%
[averts,iafaces]=qtri2atri(verts,ifaces,1);

geom_type = 2;
filename_geo_out = 'test0.a.tri';

atriwrite(filename_geo_out,averts,iafaces);
