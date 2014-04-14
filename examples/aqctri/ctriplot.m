function ctriplot(cverts,icfaces,ifrefine)
%CTRIPLOT Plot cubic triangulation.
%
%  ctriplot(cverts,icfaces,ifrefine)
%
%  Input parameters:
%
%  cverts - real(3,ncverts): array of triangulation vertices (cubic)
%  icfaces - integer(10,ncfaces): indices of triangle vertices (cubic)
%  ifrefine - if set to 0, return triangle vertices only, 
%             if set to 1, subdivide each triangle into 9, based on midpoints
%

if( nargin < 3 ), ifrefine=0; end

[verts,ifaces]=ctri2atri(cverts,icfaces,ifrefine);

%%%trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)')

nfaces=size(ifaces,2);
C=zeros(nfaces,1);
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)',C')

axis equal
axis auto
%%%alpha(.5)
