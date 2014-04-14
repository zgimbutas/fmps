function [triaquad,trianorm,triaarea,centroid,triatang1,triatang2]=ctriproc(verts,ifaces)
%CTRIPROC Process triangulations in Cart3d format (cubic).
%
%  [triaquad]=ctriproc(verts,ifaces);
%
%  [triaquad,trianorm,triaarea,centroid]=ctriproc(verts,ifaces);
%
%  [triaquad,trianorm,triaarea,centroid,triatang1,triatang2]=
%           ctriproc(verts,ifaces);
%
%  Input parameters:
%
%  verts - real(3,nverts): array of triangulation vertices
%  ifaces - integer(10,nfaces): indices of triangle vertices
%
%  Output parameters:
%
%  triaquad - real(3,10,ntri): array of triangle vertex and midpoint coordinates 
%  trianorm - real(3,nsource): triangle normals at centroids
%  triaarea - real(nsource): triangle area elements at centroids
%  centroid - real(3,nsource): triangle centroids
%  triatang1 - real(3,nsource): triangle tangents at centroids (first set)
%  triatang2 - real(3,nsource): triangle tangents at centroids (second set)
%
%  Note: the first set of tangent vectors is (\partial xyz/\partial u).
%

%
%  Construct triangle vertex array
%
nverts=size(verts,2);
nfaces=size(ifaces,2);

ntri = nfaces;
triaquad = zeros(3,10,ntri);

%for i=1:ntri	
%  triaquad(1:3,1:10,i) = verts(1:3,ifaces(1:10,i));
%end

triaquad(1:3,1:10,1:ntri) = reshape(verts(1:3,ifaces(1:10,1:ntri)),3,10,ntri);


if( nargout > 1 ),
%
%  Parametrization constants
%
%       ... setup a cubic triangle in R^3
%
%
%              2
%             . .     
%           C1   B2 
%          .       .
%         C2   M   B1
%        .           .
%       0 . A1 . A2 . 1
%
%
x0=squeeze(triaquad(1:3,1,:));
x1=squeeze(triaquad(1:3,2,:));
x2=squeeze(triaquad(1:3,3,:));
xa1=squeeze(triaquad(1:3,4,:));
xa2=squeeze(triaquad(1:3,5,:));
xb1=squeeze(triaquad(1:3,6,:));
xb2=squeeze(triaquad(1:3,7,:));
xc1=squeeze(triaquad(1:3,8,:));
xc2=squeeze(triaquad(1:3,9,:));
xm=squeeze(triaquad(1:3,10,:));

xu=-11/2*x0+x1+9*xa1-9/2*xa2;
xv=-11/2*x0+x2+9*xc2-9/2*xc1;
xuu=18*xa2+9*x0-9/2*x1-45/2*xa1;
xuv=9/4*(xa2+xc1-xb2-xb1) + 9*x0-45/4*(xc2+xa1)+27/2*xm;
xvv=18*xc1+9*x0-9/2*x2-45/2*xc2;
xuuu=9/2*(x1-x0)+27/2*(xa1-xa2);
xuuv=9*(xa1-xm)+9/2*(xc2+xb1-xa2-x0);
xuvv=9*(xc2-xm)+9/2*(xa1+xb2-xc1-x0);
xvvv=9/2*(x2-x0)+27/2*(xc2-xc1);

%
%  Triangle centroids
%
u=1/3;
v=1/3;

centroid=x0+u*xu+v*xv+u*u*xuu+2*u*v*xuv+v*v*xvv+...
   u*u*u*xuuu+3*u*u*v*xuuv+3*u*v*v*xuvv+v*v*v*xvvv;

%Test parametrization, the centroid must be equal to xm
%centroid - xm

%
%  Construct triangle normals
%
trianorm = zeros(3,ntri);
triaarea = zeros(1,ntri);

vec1 = xu+2*(u*xuu+v*xuv)+3*u*u*xuuu+6*u*v*xuuv+3*v*v*xuvv;
vec2 = xv+2*(u*xuv+v*xvv)+3*u*u*xuuv+6*u*v*xuvv+3*v*v*xvvv;

trianorm = cross(vec1,vec2);
ds = sqrt(sum(trianorm.^2,1));

trianorm = trianorm ./ repmat(ds,3,1);

%  Triangle area element at the centroid
triaarea = ds/2;
end


if( nargout > 4 ),
%
%  Construct tangent vectors
%
dt = sqrt(sum(vec1.^2,1));
triatang1 = vec1 ./ repmat(dt,3,1);
triatang2 = cross(trianorm,triatang1);
end
