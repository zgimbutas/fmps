function [a,b]=planew_ab(n)
%
% Construct multipole coefficients for incoming plane wave
% kvec = zk*[0 0 -1]; epol = [1 0 0]; 
% 
%  Input parameters:
%
%    n - number of terms in multipole expansions
%
%  Output parameters:
%
%    a,b - complex(n) - coefficients (1:n,-1) of EM-multipole expansions
%

ima=1i;
i=(1:n)';

a=sqrt(2.*i+1)/2 .* (-ima).^(i-1);
b=sqrt(2.*i+1)/2 .* (-ima).^(i);
