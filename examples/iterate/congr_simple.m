function [x,ier]=congr_simple(A,b,tol,maxit,x0)
%Conjugate gradients algorithm.
%
%  This subroutine solves a complex linear system Ax=b by means
%  of conjugate gradients algorithm.
%
%  X = CONGR_SIMPLE(A,b);
%  X = CONGR_SIMPLE(A,b,tol);
%  X = CONGR_SIMPLE(A,b,tol,maxit);
%  X = CONGR_SIMPLE(A,b,tol,maxit,x0);
%
%  Input parameters:
%
%  A - the matrix of the system (or the function handle to evaluate A(x) )
%  b - the right hand side 
%  tol - the required accuracy
%  maxit - the maximum number of iteration permitted
%  x0 - the initial guess
%
%  Output parameters:
%
%  x - the solution of the system
%  ier - error return code
%     ier=0 normal execution of the subroutine
%     ier=4 means that the maximum number iterations maxit
%           has been reached without achieving the required accuracy tol
%     ier=8 means that the errors failed to decrease before the maximum
%           number of iterations has been reached or the required accuracy
%           eps has been reached
%

[n,m]=size(b);

if( nargin < 3 ), tol = 1e-6; end
if( nargin < 4 ), maxit = min(n,20); end

ier=0;

% A is a matrix, construct function handle 
if( isnumeric(A) ), A = (@(x) A*x); end

if( nargin == 5),
x=x0;
r=b-A(x);
else
x=zeros(n,1);
r=b;
end
er0=norm(r,2);

% Randomize the initial solution vector
%x=rand(n,1);
%r=b-A(x);
%er0=norm(r,2);

p=r;
norm1 = norm(b,2); 

for i=1:maxit

% recalculate residual every 10 steps
if( mod(i,10) == 0), r=b-A(x); end;

% initial new direction guess
Ap=A(p);

%%%alpha=(r'*p) ./ (Ap'*p)
alpha=(r'*r) ./ (Ap'*p)
x=x+alpha*p;
r=r-alpha*Ap;

%%%fprintf('iter: %d, rms=norm(r,2)/sqrt(n): %17.12f\n',i,norm(r,2)/sqrt(n));
fprintf('iter: %d, rel=norm(r,2)/norm(b,2): %17.12e\n',i,norm(r,2)/norm(b,2));

%%%if( i == 1 ), norm1 = norm(r,2); end;
if( norm(r,2) < tol*norm1 ), return; end;

er1=norm(r,2);

% evaluate the quadratic form to be optimized
%q1 = x'*(b-r)/2 - b'*x;
%q1 = x'*(b-r)/2 - x'*b
%q1 = -(b+r)'*x/2
%q1 = -(b)'*x/2
%q1 = -A(x)'*x/2

% the residual stopped decreasing, abort
%if( er1 >= er0 ),
%   ier=8; return; 
%end;

beta=(er1/er0);
beta=beta^2;
p=r+beta*p;

er0=er1;

end

% no convergence, abort
ier=4;
