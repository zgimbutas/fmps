%
% script calculates moments based on the results of fmps calculations
%
% modified by Maxim Ignatenko (KSU, 06/2011)
%
function varargout = FMPS_post_mom(A,aimpole,bimpole,aompole,bompole,asmpole,bsmpole)

    eps0     = 8.8541878e-12 ;
    cmu0     = 4*pi*1e-7 ; 
    clight   = 2.99792458e8 ;
    
%     nnodes   = A.nnodes;
%     rnodes   = A.rnodes;
%     nphi     = A.nphi;
%     ntheta   = A.ntheta;
%     weights  = A.weights;

    nterms   = A.nterms;       % order of decompositions
    zk       = A.zk;           % wavenumber of incident wave
    wavelength = A.omega/2/pi; % incident wavelength [nm]
    nspheres = A.nspheres;     % number of pseusoparticles
    center   = A.center;       % positions of the centers of pseudoparticles
    radius   = A.radius;       % radii of the pseudoparticles


%
%  Evaluate the sum of electric and magnetic dipole moments 
%
    pvec_sum = zeros(3,1) + 1i*zeros(3,1);
    mvec_sum = zeros(3,1) + 1i*zeros(3,1);

    for i=1:nspheres
        [pvec0,mvec0]=empol(wavelength,zk,aompole(:,i),bompole(:,i),...
                        nterms,radius(i));	
        pvec_sum = pvec_sum + pvec0;
        mvec_sum = mvec_sum + mvec0;
    end

    pdipole_sum = pvec_sum;
    mdipole_sum = mvec_sum;


%
%  Evaluate the total electric and magnetic dipole moment 
%  with respect to x0y0z0 
%
    x0y0z0=[0 0 0]';
    [pdipole,mdipole]=empolms(wavelength,zk,nspheres,aompole,bompole,...
                        nterms,center,radius,x0y0z0);      


%
%  convert into SI units
%

%     pdipole = -pdipole*(4*pi)/Vrod/nspheres; % non-dimensionalized 
    pdipole = -pdipole*(4*pi*eps0) * 1e-27;        % SI units
    mdipole = -mdipole*(4*pi/cmu0) * 1e-27/clight; % SI units

    % to convert pdipole to non - dimensional units
    % pdipole = pdipole / ( Vrod * eps0 ) ; 


%     pdipole_sum = -pdipole_sum*(4*pi)/Vrod/nspheres; % non-dimensionalized 
    pdipole_sum = -pdipole_sum*(4*pi*eps0) * 1e-27;        % SI units
    mdipole_sum = -mdipole_sum*(4*pi/cmu0) * 1e-27/clight; % SI units

    % to convert pdipole to non - dimensional units
    % pdipole_sum = pdipole_sum / ( Vrod * eps0 ) ; 


    
    varargout{1} = pdipole;    % total dipole moments
    varargout{2} = mdipole;
    varargout{3} = pdipole_sum;    % the sum of dipole moments
    varargout{4} = mdipole_sum;

end
