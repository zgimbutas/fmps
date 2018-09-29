function [frames]=em3orient(rot,nspheres)
%EM3ORIENT Convert Euler angles into orientation frames.
%
%  FRAMES=EM3ORIENT(SPHERE_ROT,NSPHERES);
%
%  Input parameters:
%
%    sphere_rot - real(3,nspheres): Euler angles for a collection of spheres
%    nspheres - the number of spheres
%
%  Output parameters:
%
%    frames - real(3,3,nspheres): rotation matrices (frames) for all spheres
%

frames = zeros(9,nspheres);

mex_id_ = 'em3orient(i double[], i int[x], io double[])';
[frames] = fmps(mex_id_, rot, nspheres, frames, 1);

frames = reshape(frames,3,3,nspheres);



