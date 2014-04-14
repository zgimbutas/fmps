function [sphere_rot]=euler_angles(itype,nspheres)
%EULER_ANGLES Construct Euler rotation angles standard orientation types.
% 
%  [sphere_rot]=euler_angles(itype,nspheres);
%
%  Input parameters:
% 
%    itype - rotation type
%      itype=1: standard orientation, z-axis
%      itype=2: x-axis orientation (z -> x, x -> -z, y -> y)
%      itype=3: y-axis orientation (z -> -y, x -> -z, y -> x)
%      itype=4: x-axis orientation (z -> x, x -> -y, y -> -z)
%      itype=5: y-axis orientation (z -> -y, x -> -x, y -> -z)
%      itype=6: random orientation
%      itype=7: pseudo-random orientation with twists
%      itype=8: z-rotation (x -> -y, y -> x, z -> z)
%
%    nspheres - the number of spheres
%
%  Output parameters:
% 
%    sphere_rot - real(3,nspheres): Euler angles for a collection of spheres
%

sphere_rot=zeros(3,nspheres);

%       
%       ... Euler rotation angles for all inclusions
%
for kk=1:nspheres
%
%       ... standard orientation, z-axis
%
if( itype == 1 )
   sphere_rot(1,kk)=0;
   sphere_rot(2,kk)=0;
   sphere_rot(3,kk)=0;
end
%
%       ... x-axis orientation (z -> x, x -> -z, y -> y)
%
if( itype == 2 )
   sphere_rot(1,kk)=0;
   sphere_rot(2,kk)=pi/2;
   sphere_rot(3,kk)=0;
end
%
%       ... y-axis orientation (z -> -y, x -> -z, y -> x)
%
if( itype == 3 )
   sphere_rot(1,kk)=pi/2;
   sphere_rot(2,kk)=pi/2;
   sphere_rot(3,kk)=0;
end
%
%        ... x-axis orientation (z -> x, x -> -y, y -> -z)
%
if( itype == 4 )
   sphere_rot(1,kk)=0;
   sphere_rot(2,kk)=pi/2;
   sphere_rot(3,kk)=pi/2;
end
%
%        ... y-axis orientation (z -> -y, x -> -x, y -> -z)
%
if( itype == 5 )
   sphere_rot(1,kk)=pi/2;
   sphere_rot(2,kk)=pi/2;
   sphere_rot(3,kk)=pi/2;
end
%
%       ... random orientation
%
if( itype == 6 )
   sphere_rot(1,kk)=pi*2*rand;
   sphere_rot(2,kk)=pi*rand;
   sphere_rot(3,kk)=pi*2*rand;
end
%
%       ... pseudo-random orientation with twists
%
if( itype == 7 )
   sphere_rot(1,kk)=pi*2*(kk)/nspheres;
   sphere_rot(2,kk)=pi*(kk)/nspheres;
   sphere_rot(3,kk)=pi*2*(kk)/nspheres;
end
%
%       ... z-rotation (x -> -y, y -> x, z -> z)
%
if( itype == 8 )
   sphere_rot(1,kk)=pi/2;
   sphere_rot(2,kk)=0;
   sphere_rot(3,kk)=0;
end
%
end

