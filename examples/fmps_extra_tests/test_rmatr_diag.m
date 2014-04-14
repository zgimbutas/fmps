
%
%  Construct reflection coefficients for arbitrary dielectric sphere
%
test_die_arb

%
%  Construct the diagonal of reflection matrix
%  Not the most efficient way to do that but it works
%

%x = (1:(nterms+1)).^2-1
%y = diff(x)
y = 2*(1:nterms) + 1

%
%  Set reflection coefficient for zero'th mode to zero
%
raa_diag = 0;
rbb_diag = 0;
for i=1:nterms
  raa_diag = [raa_diag; repmat(ra(i),y(i),1)];
  rbb_diag = [rbb_diag; repmat(rb(i),y(i),1)];
end

raa_diag
rbb_diag

size(raa_diag);
size(rbb_diag);
