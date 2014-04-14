function [verts,ifaces]=qtri2atri(qverts,iqfaces,ifrefine)
%QTRI2ATRI Convert quadratic triangulation to flat triangulation.
%
%  [verts,ifaces]=qtri2atri(qverts,iqfaces)
%  [verts,ifaces]=qtri2atri(qverts,iqfaces,ifrefine)
%
%  Input parameters:
%
%  qverts - real(3,nqverts): array of triangulation vertices (quadratic)
%  iqfaces - integer(6,nqfaces): indices of triangle vertices (quadratic)
%  ifrefine - if set to 0, return triangle vertices only, 
%             if set to 1, subdivide each triangle into 4, based on midpoints
%
%  Output parameters:
%
%  verts - real(3,nverts): array of triangulation vertices (flat)
%  ifaces - integer(3,nfaces): indices of triangle vertices (flat)
%

%
%  Convert quadratic triangulation to flat triangulation,
%  while preserving orientation 
%

if( nargin < 3 ), ifrefine = 0; end


if( ifrefine == 0 ),

nqverts=size(qverts,2);
nqfaces=size(iqfaces,2);

nverts = nqverts;
nfaces = nqfaces;

verts = qverts;

ifaces = [iqfaces(1,:); iqfaces(2,:); iqfaces(3,:)];

end


if( ifrefine == 1 ),

nqverts=size(qverts,2);
nqfaces=size(iqfaces,2);

nverts = nqverts;
nfaces = 4*nqfaces;

verts = qverts;

ifaces = [[iqfaces(1,:); iqfaces(4,:); iqfaces(6,:)], ...
          [iqfaces(2,:); iqfaces(5,:); iqfaces(4,:)], ...
          [iqfaces(3,:); iqfaces(6,:); iqfaces(5,:)], ...
          [iqfaces(4,:); iqfaces(5,:); iqfaces(6,:)]];

end

