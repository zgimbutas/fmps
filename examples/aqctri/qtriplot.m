function qtriplot(qverts,iqfaces,ifrefine)
%QTRIPLOT Plot quadratic triangulation.
%
%  qtriplot(qverts,iqfaces,ifrefine)
%
%  Input parameters:
%
%  qverts - real(3,nqverts): array of triangulation vertices (quadratic)
%  iqfaces - integer(6,nqfaces): indices of triangle vertices (quadratic)
%  ifrefine - if set to 0, return triangle vertices only, 
%             if set to 1, subdivide each triangle into 4, based on midpoints
%

if( nargin < 3 ), ifrefine=0; end

[verts,ifaces]=qtri2atri(qverts,iqfaces,ifrefine);

%%%trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)')

nfaces=size(ifaces,2);
C=zeros(nfaces,1);
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)',C')

axis equal
axis auto
%%%alpha(.5)
