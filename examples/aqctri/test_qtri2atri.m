
%
%  Retrieve quadratic triangulation
%

geom_type = 3;
%%%filename_geo = 'test1.q.tri';  
%%%filename_geo = 'test2.q.tri';
filename_geo = 'sphere80.q.tri';

[verts,ifaces,nverts,nfaces] = qtriread(filename_geo);
nverts,nfaces

%
%  Construct triangle vertex, normal, area, and centroid arrays
%

ntri = nfaces;
[triangles,trianorm,triaarea,source]=qtriproc(verts,ifaces);


%
%  Convert quadratic triangulation in flat one, and store
%
[averts,iafaces]=qtri2atri(verts,ifaces,1);

geom_type = 2;
filename_geo_out = 'test1.a.tri';

atriwrite(filename_geo_out,averts,iafaces);
