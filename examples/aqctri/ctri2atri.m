function [verts,ifaces]=ctri2atri(cverts,icfaces,ifrefine)
%CTRI2ATRI Convert cubic triangulation to flat triangulation.
%
%  [verts,ifaces]=qtri2atri(cverts,icfaces)
%  [verts,ifaces]=qtri2atri(cverts,icfaces,ifrefine)
%
%  Input parameters:
%
%  cverts - real(3,ncverts): array of triangulation vertices (cubic)
%  icfaces - integer(10,ncfaces): indices of triangle vertices (cubic)
%  ifrefine - if set to 0, return triangle vertices only, 
%             if set to 1, subdivide each triangle into 9, based on midpoints
%
%  Output parameters:
%
%  verts - real(3,nverts): array of triangulation vertices (flat)
%  ifaces - integer(3,nfaces): indices of triangle vertices (flat)
%

%
%  Convert cubic triangulation to flat triangulation,
%  while preserving orientation 
%

if( nargin < 3 ), ifrefine = 0; end


if( ifrefine == 0 ),

ncverts=size(cverts,2);
ncfaces=size(icfaces,2);

nverts = ncverts;
nfaces = ncfaces;

verts = cverts;

ifaces = [icfaces(1,:); icfaces(2,:); icfaces(3,:)];

end


if( ifrefine == 1 ),

ncverts=size(cverts,2);
ncfaces=size(icfaces,2);

nverts = ncverts;
nfaces = 10*ncfaces;

verts = cverts;

ifaces = [[icfaces(1,:); icfaces(4,:); icfaces(9,:)], ...
          [icfaces(4,:); icfaces(5,:); icfaces(10,:)], ...
          [icfaces(5,:); icfaces(2,:); icfaces(6,:)], ...
          [icfaces(6,:); icfaces(7,:); icfaces(10,:)] ...
          [icfaces(7,:); icfaces(3,:); icfaces(8,:)], ...
          [icfaces(8,:); icfaces(9,:); icfaces(10,:)], ...
          [icfaces(9,:); icfaces(4,:); icfaces(10,:)], ...
          [icfaces(5,:); icfaces(6,:); icfaces(10,:)], ...
          [icfaces(7,:); icfaces(8,:); icfaces(10,:)]];

end

