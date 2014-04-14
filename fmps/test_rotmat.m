%
% Test Euler angle routines
%

% arbitrary rotation
sphere_rot=[1 -2 3];
sphere_rot=(rand(1,3)-.5)*pi;

%      itype=2: x-axis orientation (z -> x, x -> -z, y -> y)
%      itype=3: y-axis orientation (z -> -y, x -> -z, y -> x)
%      itype=4: x-axis orientation (z -> x, x -> -y, y -> -z)

%sphere_rot=[0 pi/2 0];  % itype=2
%sphere_rot=[pi/2 pi/2 0];  % itype=3
%sphere_rot=[0 pi/2 pi/2];  % itype=4

sphere_rot


%%%rotmat = em3orient(sphere_rot(:,i),1);

%%% Reconstruct the frame rotation matrix from Euler angles
rotmat1 = em3orient(sphere_rot,1)

%%% The same thing in Matlab
gamma = sphere_rot(3);
beta  = sphere_rot(2);
alpha = sphere_rot(1);

rot1 = [cos(gamma) sin(gamma) 0;
       -sin(gamma) cos(gamma) 0;
        0          0          1];

rot2 = [cos(beta) 0 sin(beta);
        0         1        0 ;
       -sin(beta) 0 cos(beta)];

rot3 = [cos(alpha) sin(alpha) 0;
       -sin(alpha) cos(alpha) 0;
        0          0          1];

rotmat2 = rot3 * rot2 * rot1

%%% Check the error
error=norm(rotmat1-rotmat2,2)

if( 1 == 2 ),
%%% Print the rotated frame 0x'y'z' coordinates
xprime=rotmat2*[1 0 0]'
yprime=rotmat2*[0 1 0]'
zprime=rotmat2*[0 0 1]'
end



% Reconstruct the Euler angles
euler=em3euler(rotmat2)
%%% Check the error in reconstructing the Euler angles
error=norm(euler-sphere_rot,2)

% Test the alternative form of Euler angles
euler_alt=em3euler_alt(euler,1)
%%% Check the error in reconstructing the Euler angles
error=norm(euler_alt-sphere_rot,2)

