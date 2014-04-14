function [raa_diag,rbb_diag]=rmatr_diag(nterms,ra,rb)
%
%  Construct the diagonally scaled reflection matrix.
%

y = 2*(1:nterms)+1;

%
%  Set reflection coefficient for zero'th mode to zero
%
raa_diag = 0;
rbb_diag = 0;
for i=1:nterms
  raa_diag = [raa_diag; repmat(ra(i),y(i),1)];
  rbb_diag = [rbb_diag; repmat(rb(i),y(i),1)];
end

