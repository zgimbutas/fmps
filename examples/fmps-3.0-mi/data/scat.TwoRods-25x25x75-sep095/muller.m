%
%  EM scattering from multiple dielectric and PEC spheres 
%  EM scattering from multiple dielectric inclusions
%
%  gold_jc_data contains a table of the complex refractive index for gold
%  at various frequencies (Johnson and Christy, Phys Rev B, 1972).
%
%  Column 1 is energy (eV), Column 2 is wavelength in nm,
%  Column 3 is the real part of the refractive index (n).
%  Column 4 is the imaginary part of the refractive index (k).
%
%
%  Step 1: Read in material parameters for desired frequency/frequencies
%

% updated/edited by G.P.R.  

tic

% set parameters
gold_jc_data  ; % loads the refraction index for Au at various wavelengths
freqMin   = 40; % Lowest index (longest wavelength) of gold data to be used
freqMax   = 70; % Highest index of gold data to be used
% radius    = 40; % radius of the enclosing sphere
radius    = 0.5*sqrt((25+95)^2+75^2) + 1; % radius of the enclosing sphere
nterms    = 5 ; % order of some expansion
geom_type = 2 ; % 2 = flat triangulation, 3 = quadratic triangulation

GMREStol  = 1e-06 ; % solver tolerance
GMRESinit = 20    ; % number of iterations before forced restart...

% filename_geo = 'TwoRods-20x20x50-sep30.a.tri' ; % filename with triangle list
filename_geo = 'TwoRods-25x25x75-sep95.a.tri';

%
%
%  Step 2:  Set frequency scan (using data read in from gold_jc_data).
%           Get wavelength (in nm).
%           For each wavelength, set complex refractive index (re_n,im_n).
%           Set geometry type. 
%	            For a flat triangulation, geom_type = 2
%               For a quadrati c triangulation, geom_type = 3
%               See release notes appendix for more information.
%           Get discretization by setting filename_geo to the file containing
%               the discretization.
%               NOTE: The structure in filename_geo is assumed to 
%                     centered at the origin.
%           Stretch or shift location of inclusions.
%               scale_geo=[a b c] stretches the x,y,z coordinates 
%               as read from filename_geo.
%               shift_geo=[a b c] moves the structure by the vector [a,b,c].
%           Set radius of enclosing sphere (in nm).
%           Set order of expansion 'nterms'.
%           Set solver_type = 2 for GMRES iterative solution to integral equation.
%           Set solver tolerance eps. 
%           Set GMRES parameter numit = 40 (restart parameter in GMRES). Do not modify
%               without some GMRES expertise.
%           Set output files 'filename_out' where scattering matrices will be written.
%
for freq = freqMin:freqMax

  wavelength=gold_jc(freq,2);
  re_n=gold_jc(freq,3)    ; % real part of reflective index at current frequency
  im_n=gold_jc(freq,4)    ; % imag  ""      ""        ""        ""      ""

  scale_geo=[1, 1, 1]     ; % scaling of triangles in x,y,z directions
  shift_geo=[0, 0, 0]     ; % shift of triangles in the x,y,z directions

  solver_type = 2         ; % 2 = GMRES 3 = BGCG-Stab
  eps         = GMREStol  ; % solver tolerance
  numit       = GMRESinit ; % number of iterations before restart (I think)

  filename_out = ['scat.freq.' num2str(freq)];  % output file

%
%   A number of parameters set above are needed by the underlying Fortran-based
%   solvers. They are written to the file 'config.freq' here.
%

  config = ['config.freq.' num2str(freq)];

  fid = fopen(config,'w');

  fprintf(fid,'%d\n',geom_type);
  fprintf(fid,'%s\n',filename_geo);
  fprintf(fid,'%g %g %g %g %g %g\n', scale_geo(1:3)', shift_geo(1:3)');
  fprintf(fid,'%g\n', radius);
  fprintf(fid,'%g %g %g\n', wavelength, re_n, im_n);
  fprintf(fid,'%d\n', nterms);
  fprintf(fid,'%d %g %d\n', solver_type,eps,numit);
  fprintf(fid,'%s\n',filename_out);

  fclose(fid);

%
%   Now call the solver interface code (using the appropriate Operating System).
%   int2_w32 for Windows, 32 bit
%   int2_w32_omp for Windows on multicore 32 bit machine with OpenMP
%   int2_w64 for Windows, 64 bit
%   int2_lnx64_omp for Linux, 64 bit with OpenMP
%   int2_macosi for Mac OS, Intel 32 bit 
%

% system(['int2_w32.exe ' config]);
% system(['int2_w32_omp.exe ' config]);
% system(['int2_w64.exe ' config]);
system(['./int2_lnx64_omp ' config]);
% system(['./int2_macosi ' config]);

end

t=toc;
fprintf('running time is %f sec\n',t);