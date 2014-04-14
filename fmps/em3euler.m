function euler=em3euler(rotmat,ifunwrap)
%EM3EULER Reconstruct Euler angles from a rotation matrix
%
% Given a rotation matrix, reconstruct the Euler angles.
% This is an inverse of em3orient.
%
%  EULER=EM3EULER(ROTMAT);
%  EULER=EM3EULER(ROTMAT,IFUNWRAP);
%
%  Input parameters:
%
%    rotmat - real(3,3): rotation matrix (frame) 
%    ifunwrap - if set to 1, unwrap the phase
%
%  Output parameters:
%
%    euler - real(3): Euler angles 
%
%
% Note 1: Euler angles are not unique, complicated phase wrapping may occur.
%
% Note 2: Our definition of Euler angles has a negative sign for beta angle
% for compatibility with rotviarecur3 and rotproj routines.
%

if( nargin < 2 ), ifunwrap = 0; end

% Recover the last two Euler angles, 
x=rotmat(3,1);
y=rotmat(3,2);
z=rotmat(3,3);

if( abs(x) == 0 && abs(y) == 0 ),
  phi = 0;
else
  phi = atan2(y,x);
end

theta=acos(z);

beta=-theta;
gamma=phi;

% We have recovered two Euler angles, 
% rotate into the new coordinate system to recover
% the remaining rotation about z axis 

rot1 = [cos(gamma) sin(gamma) 0;
       -sin(gamma) cos(gamma) 0;
        0          0          1];

rot2 = [cos(beta) 0 sin(beta);
        0         1        0 ;
       -sin(beta) 0 cos(beta)];

%%% recover rot3 by using  rotmat = rot3 * rot2 * rot1;

rotmat_alpha=rotmat*rot1'*rot2';

x=rotmat_alpha(1,1);
y=rotmat_alpha(1,2);

if( abs(x) == 0 && abs(y) == 0 ),
  phi_alpha = 0;
else
  phi_alpha = atan2(y,x);
end

alpha=phi_alpha;

rot3 = [cos(alpha) sin(alpha) 0;
       -sin(alpha) cos(alpha) 0;
        0          0          1];


euler=[alpha -theta phi];
%%%euler=[alpha beta gamma];


if( ifunwrap == 1 ),

% yet another way to standartize Euler angles, always return beta>0
%if( -theta < 0 ),
%  euler=[alpha+pi theta phi-pi]
%end

% unwrap the phase
if( euler(1) > +pi ), euler(1)=euler(1)-2*pi; end
if( euler(1) < -pi ), euler(1)=euler(1)+2*pi; end
if( euler(3) > +pi ), euler(3)=euler(3)-2*pi; end
if( euler(3) < -pi ), euler(3)=euler(3)+2*pi; end
end


% optionally, check the error 
%%%error_euler_rot=norm(rotmat - rot3*rot2*rot1,2)
