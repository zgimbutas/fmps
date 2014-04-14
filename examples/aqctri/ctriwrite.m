function [nverts,nfaces]=ctriwrite(filename,verts,ifaces)
%CTRIWRITE Store triangulations in Cart3d to a file.  (cubic)
%
%  [nverts,nfaces]=ctriwrite(filename,verts,ifaces);
%
%  Input parameters:
%
%  filename - output file name.
%  verts - real(3,nverts): array of triangulation vertices
%  ifaces - integer(10,nfaces): indices of triangle vertices
%
%  Output parameters:
%
%  nverts - number of triangulation vertices
%  nfaces - number of triangle vertex indices 
%

%
%  Store cubic triangulation
%

geom_type = 4;

fid = fopen(filename,'w');

nverts=size(verts,2);
nfaces=size(ifaces,2);

fprintf(fid,'%d %d\n',nverts,nfaces);

fprintf(fid,'%20.15e %20.15e %20.15e\n',verts);
fprintf(fid,'%d %d %d %d %d %d %d %d %d %d\n',ifaces);

fclose(fid);

%%%nverts,nfaces

