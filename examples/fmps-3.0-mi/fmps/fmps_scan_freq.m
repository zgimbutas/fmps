%
%  Sample EM multi-sphere scattering code
%  scan over different frequencies but the same pseudoparticle
%
% modified by Maxim Ignatenko (KSU, 06/2011)
%
clear variables

% folder with utilites (functions, mex-files ...)
UtilFolderName = [pwd filesep 'utils'];
DataFolderName = ['..' filesep 'data'];
InfoFileName   = 'info.txt';
addpath(UtilFolderName); % add folder to MatLab path

% frequently changed params
NumSpheres = 2;
% list of folders with results (these folders are in data folder)
if NumSpheres == 1
    NumRods        = 2 ;    % number of rods in a pseudoparticle
    FileFolderName = 'scat.TwoRods-25x25x75-sep095';
else
    NumRods        = 1 ;    % number of rods in a pseudoparticle
    FileFolderName = 'scat.OneRod-25x25x75'; % scattering data location
end;


% precalculated separation between rods
SepDist       = 95;

FreqMin       = 40; % lowest index (longest wavelength) of gold data to be tested
FreqMax       = 70; % highest index of gold data to be tested

PropDir       = [1 0 0]; % direction of propagation of incident field
Polarization  = [0 0 1];  % polarization of incident E-field

% NumRods       = 2 ;    % number of rods in a pseudoparticle
RodRadius     = 12.5 ; % radius of rod for volume calculation
AspectRatio   = 3 ;    % aspect ratio of rod for volume calculation

% generate vector with numbers of frequencies
Freq          = FreqMin:FreqMax;

fprintf('\n')
fprintf('*********************************************\n');
fprintf('START OF CALCULATIONS\n');
fprintf('*********************************************\n');

Pdipole_1d = zeros(3, length(Freq));
Mdipole_1d = zeros(3, length(Freq));
Qabs_1d    = zeros(1, length(Freq));
Qsca_1d    = zeros(1, length(Freq));

    
    % automatically take some parameters
    InfoFile      = [DataFolderName filesep ...
                     FileFolderName filesep InfoFileName];

    fid           = fopen(InfoFile);
    data          = textscan(fid,'%s');
    fclose(fid);

    tmp = data{1}(9);  SphereRadius = sscanf(tmp{1},'%d');% sphere radii
    tmp = data{1}(13); nterms       = sscanf(tmp{1},'%d');% order of decompositions
    
%     nterms = 1;

    %  Set sphere type in array sphere_type 
    %
    %      For type=1, the sphere is a perfect conductor.
    %      For type=2, the sphere is a dielectric sphere.
    %      For type=3, use scattering matrix obtained from Muller solver.
    %  Set permeability and permittivity for the inclusions if itype = 2.
    %      (ignored for itype = 1,3)
    %
    SphereType = 3;
    SphereEpr  = 1; % used only for SphereType = 2
    SphereMur  = 1; % used only for SphereType = 2
    
if NumSpheres == 1
    % init geometry with one sphere at the origin
    A          = FMPS_init_single(nterms, SphereRadius, SphereType, ...
        SphereEpr, SphereMur);
else
    % init geometry with two spheres (symmetrical relative origin)
    A          = FMPS_init_two(nterms, SphereRadius, SphereType, ...
        SphereEpr, SphereMur, SepDist*[1 0 0]);
end;

    % read optical parameters
    gold_jc_data;

    % set FMM precision parameter
    % A.iprec = 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %  STEP 4. Run frequency scan (using gold_jc_data points)
    %

    for ifreq = 1:length(Freq)

        %
        % prepare to run FMPS
        %

        % take wavelength [nm]
        wavelength = gold_jc(Freq(ifreq),2) ;
        
        % copy precalculated muller data to 'fort.19'
        scatfile   = [DataFolderName filesep FileFolderName ...
              filesep 'scat.freq.' num2str(Freq(ifreq))] ; 
        copyfile(scatfile,'fort.19');

        %
        % run FMPS
        %
        [aimpole,bimpole,aompole,bompole,asmpole,bsmpole,A] = ...
            FMPS_run(A, wavelength, PropDir, Polarization) ;

        %
        % postprocess FMPS
        %
        [qabs,qsca,qext] = ...
            FMPS_post_Q(A, aimpole,bimpole,aompole,bompole,asmpole,bsmpole) ;
        [pdipole,mdipole] = ...
            FMPS_post_mom(A, aimpole,bimpole,aompole,bompole,asmpole,bsmpole) ;

        Qabs_1d(ifreq)      = qabs;   % [W]
        Qsca_1d(ifreq)      = qsca;   % [W]
        Pdipole_1d(:,ifreq) = pdipole;% in units SI
        Mdipole_1d(:,ifreq) = mdipole;% in units SI

    end;

fprintf('*********************************************\n');
fprintf('END OF CALCULATIONS\n');
fprintf('*********************************************\n');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   PLOT RESULTS
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % if dimensionless pdipole is desirable
% eps0     = 8.8541878e-12 ;
% cmu0     = 4*pi*1e-7 ; 
% clight   = 2.99792458e8 ;
% % volume of rods in psuedo-particle (for the estimation of moments)
% Vrod     = NumSpheres*NumRods * pi * RodRadius^3 * ...
%     ( 4 / 3  + 2 * (AspectRatio - 1));
% Vrod     = Vrod*1e-27; % nm^3 -> m^3
% % convert to dimensioneless 
% Pdipole_1d = Pdipole_1d/eps0/Vrod;

%
%  plot electric polarization vector
%
figure;
fplot = Pdipole_1d;
x     = gold_jc(Freq,2);

subplot(3,1,1)
hold on
plot( x, real(fplot(1,:)), 'b', 'LineWidth', 2)
plot( x, imag(fplot(1,:)), 'r', 'LineWidth', 2)
plot( x, zeros(size(x)), 'k--');
hold off
xlabel('wavelength [nm]')
ylabel('x-direction')
title(['electric polarization vector (SI), separation = ', ...
    num2str(SepDist), 'nm'])
legend( 'Real', 'Imag','Location','SouthEast' ) ; 

subplot(3,1,2)
hold on
plot( x, real(fplot(2,:)), 'b', 'LineWidth', 2)
plot( x, imag(fplot(2,:)), 'r', 'LineWidth', 2)
plot( x, zeros(size(x)), 'k--');
hold off
xlabel('wavelength [nm]')
ylabel('y-direction')
legend( 'Real', 'Imag','Location','SouthEast' ) ; 

subplot(3,1,3)
hold on
plot( x, real(fplot(3,:)), 'b', 'LineWidth', 2)
plot( x, imag(fplot(3,:)), 'r', 'LineWidth', 2)
plot( x, zeros(size(x)), 'k--');
hold off
xlabel('wavelength [nm]')
ylabel('z-direction')
legend( 'Real', 'Imag','Location','SouthEast' ) ; 

%
%  plot magnetic polarization vector
%
figure;
fplot = Mdipole_1d;
x     = gold_jc(Freq,2);

subplot(3,1,1)
hold on
plot( x, real(fplot(1,:)), 'b', 'LineWidth', 2)
plot( x, imag(fplot(1,:)), 'r', 'LineWidth', 2)
plot( x, zeros(size(x)), 'k--');
hold off
xlabel('wavelength [nm]')
ylabel('x-direction')
title(['magnetic polarization vector (SI), wavelength = ', ...
    num2str(wavelength), 'nm'])
legend( 'Real', 'Imag','Location','SouthEast' ) ; 

subplot(3,1,2)
hold on
plot( x, real(fplot(2,:)), 'b', 'LineWidth', 2)
plot( x, imag(fplot(2,:)), 'r', 'LineWidth', 2)
plot( x, zeros(size(x)), 'k--');
hold off
xlabel('wavelength [nm]')
ylabel('y-direction')
legend( 'Real', 'Imag','Location','SouthEast' ) ; 

subplot(3,1,3)
hold on
plot( x, real(fplot(3,:)), 'b', 'LineWidth', 2)
plot( x, imag(fplot(3,:)), 'r', 'LineWidth', 2)
plot( x, zeros(size(x)), 'k--');
hold off
xlabel('wavelength [nm]')
ylabel('z-direction')
legend( 'Real', 'Imag','Location','SouthEast' ) ; 

figure
fplot = Pdipole_1d;
x     = gold_jc(Freq,2);
hold on
plot( x, real(fplot(3,:)), 'b', 'LineWidth', 2)
plot( x, imag(fplot(3,:)), 'r', 'LineWidth', 2)
plot( x, zeros(size(x)), 'k--');
hold off
xlabel('wavelength [nm]')
ylabel('electric polarization vector, z-direction')
title(['separation = ', num2str(SepDist), 'nm'])
legend( 'Real', 'Imag','Location','SouthEast' ) ; 
% saveas(gcf,'ElecPolVector.fig') 

figure
fplot = Mdipole_1d;
x     = gold_jc(Freq,2);
hold on
plot( x, real(fplot(2,:)), 'b', 'LineWidth', 2)
plot( x, imag(fplot(2,:)), 'r', 'LineWidth', 2)
plot( x, zeros(size(x)), 'k--');
hold off
xlabel('wavelength [nm]')
ylabel('magnetic polarization vector, y-direction')
title(['separation = ', num2str(SepDist), 'nm'])
legend( 'Real', 'Imag','Location','SouthEast' ) ; 
% saveas(gcf,'MagPolVector.fig') 

%
%  plot absorption and scattering x-sections
%
figure;
x     = gold_jc(Freq,2);

hold on
plot( x, -real(Qabs_1d), 'b', 'LineWidth', 2)
plot( x,  real(Qsca_1d), 'r', 'LineWidth', 2)
hold off
xlabel('wavelength [nm]')
ylabel('absorption and scattering cross sections')
title(['separation = ', num2str(SepDist), 'nm'])
legend( 'Abs', 'Scat','Location','SouthEast' ) ; 
grid on
% saveas(gcf,'AbsorptionCrossSections.fig')

%
%  POSTPROCESSING - part 2
%
%  evaluate E and H fields at a target grid
%   1 wavelengths below the x-y plane
%

% % % ngrid = 64;
% % % 
% % % x = linspace(-2100,2100,ngrid);
% % % y = linspace(-2100,2100,ngrid);
% % % 
% % % %evecs = zeros(3,ngrid,ngrid);
% % % %hvecs = zeros(3,ngrid,ngrid);
% % % 
% % % ntargets = ngrid*ngrid;
% % % targets = zeros(3,ngrid,ngrid);
% % % for k=1:ngrid
% % %     for j=1:ngrid
% % %         targets(:,j,k) = [x(j), y(k), -wavelength]';
% % %     end
% % % end
% % % 
% % % % scattered field
% % % [evecs,hvecs] = em3d_mpole_targeval(A.nspheres,A.nterms,A.ncoefs,...
% % %     A.omega,A.eps0,A.cmu0,A.center,A.radius,...
% % %     aompole,bompole,ntargets,targets);
% % % 
% % % % incoming field
% % % kvec = zk*[0 0 -1];
% % % epol = [1 0 0];
% % % [evec0,hvec0] = em3d_planearb_targeval(kvec,epol,ntargets,targets);
% % % 
% % % % total field
% % % evecs = evecs+evec0;
% % % hvecs = hvecs+hvec0;
% % % 
% % % evecs = reshape(evecs,3,ngrid,ngrid);
% % % hvecs = reshape(hvecs,3,ngrid,ngrid);


% figure(11)
% imagesc(x,y,reshape(real(evecs(1,:,:)),ngrid,ngrid))
% colorbar
% title('Total field: Re(E_x)')
% saveas(gcf,'TotalFieldRealEx.fig') ;
% 
% 
% figure(12)
% imagesc(x,y,reshape(imag(evecs(1,:,:)),ngrid,ngrid))
% colorbar
% title('Total field: Im(E_x)')
% saveas(gcf,'TotalFieldImagEx.fig') ; 
% 
% figure(13)
% imagesc(x,y,reshape(angle(evecs(1,:,:)),ngrid,ngrid))
% colorbar
% title('Total field: phase(E_x)')
% saveas(gcf,'TotalFieldPhaseEx.fig') 

% cause I forget to save things all the time
%save SCANvariables.mat

% remove pathes from MatLab path
rmpath(UtilFolderName);

