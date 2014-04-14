function [verts,ifaces,nverts,nfaces]=ctriread(filename)
%CTRIREAD Retrieve triangulations in Cart3d from a file.  (cubic)
%
%  [verts,ifaces,nverts,nfaces]=ctriread(filename);
%
%  Input parameters:
%
%  filename - input file name.
%
%  Output parameters:
%
%  verts - real(3,nverts): array of triangulation vertices
%  ifaces - integer(10,nfaces): indices of triangle vertices
%

%
%  Retrieve cubic triangulation
%

geom_type = 4;

fid = fopen(filename,'r');

nverts=0;
nfaces=0;
[nverts] = fscanf(fid,'%d',1);
[nfaces] = fscanf(fid,'%d',1);

verts=zeros(3,nverts);
ifaces=zeros(3,nfaces);

[verts] = fscanf(fid,'%f',[3,nverts]);
[ifaces] = fscanf(fid,'%d',[10,nfaces]);

fclose(fid);

%%%nverts,nfaces

