function fmps_prini(unit1,unit2)
%FMPS_PRINI Initialize internal printing/debugging routines.
%
% Calling FMPS_PRINI(6,13) causes printing to screen and file fort.13.     
%

if (nargin == 1 )
unit2=0;
end

mex_id_ = 'prini(i int[x], i int[x])';
fmps_r2012a(mex_id_, unit1, unit2, 1, 1);

