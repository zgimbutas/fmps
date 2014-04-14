function euler_alt=em3euler_alt(euler,ifunwrap)
%EM3EULER_ALT Alternative representation of Euler angles.
%
%  Return a second form of Euler angles with flipped beta angle, 
%  and reflected alpha and gamma.
%
%  Note, that Euler angles are not unique, complicated phase wrapping may occur.
%  The alternative representation generates the same rotation matrix.
%
%  EULER_ALT=EM3EULER_ALT(EULER);
%  EULER_ALT=EM3EULER_ALT(EULER,IFUNWRAP);
%
%  Input parameters:
%
%    euler - real(3): Input Euler angles 
%    ifunwrap - if set to 1, unwrap the phase
%
%  Output parameters:
%
%    euler_alt - real(3): Alternative Euler angles 
%
%

euler_alt=[euler(1)+pi -euler(2) euler(3)-pi];

if( ifunwrap == 1 ),
% unwrap the phase
if( euler_alt(1) > +pi ), euler_alt(1)=euler_alt(1)-2*pi; end
if( euler_alt(1) < -pi ), euler_alt(1)=euler_alt(1)+2*pi; end
if( euler_alt(3) > +pi ), euler_alt(3)=euler_alt(3)-2*pi; end
if( euler_alt(3) < -pi ), euler_alt(3)=euler_alt(3)+2*pi; end
end

