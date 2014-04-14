%
% script initializes fmps software
%
% this is a simple geometry:
% a single psuedoparticle is located at the center point [0, 0, 0],
% the geometry of the pseudoparticle is kept at original state
%

% modified by Maxim Ignatenko (KSU, 06/2011)

function varargout = FMPS_init_two(nterms, SphereRadius, SphereType, ...
    SphereEpr, SphereMur, SepDist)

%
%  Sample EM multi-sphere scattering code
%  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  STEP 1
%
%  First, set solver parameters.
%
%  Set number of terms in multipole expansions 
%      (must be less than or equal to the order of the scattering 
%       matrix generated by the Muller solver). 
%
%  Set itype  (internal parameter, leave set to 1).
%
%  Set nquad  (leave set at nterms    - used internally by Fortran code)
%  Set nphi   (leave set at 2*nquad+1 - used internally by Fortran code)
%  Set ntheta (leave set at nquad+1   - used internally by Fortran code)
%  Set A as below.

%    % solver params
    itype      = 1;
    nquad      = nterms;
    nphi       = 2*nquad+1;
    ntheta     = nquad+1;

    A.nterms   = nterms;
    A.nquad    = nquad;
    A.nphi     = nphi;
    A.ntheta   = ntheta;

    [A.rnodes,A.weights,A.nnodes] = e3fgrid(itype,nquad,nphi,ntheta);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  STEP 2
%
%  Set geometry information 
%
%  Set number of spheres 
%  Set sphere centers in array center
%  Set sphere radii in array radius
%  Set sphere type in array sphere_type 
%      For type=1, the sphere is a perfect conductor.
%      For type=2, the sphere is a dielectric sphere.
%      For type=3, use scattering matrix obtained from Muller solver.
%  Set permeability and permittivity for the inclusions if itype = 2.
%      (ignored for itype = 1,3)
%  Set sphere_rot to change orientation of inclusions if itype = 3.
%  sphere_rot(1:3,i) = [0 0 0] means leave orientation as that from which the 
%                      muller solver computed the scattering matrix 
%  sphere_rot(1:3,i) defines the 3 Euler angles that determine an arbitrary orientation.
%  The file rotangles.m contains some typical examples of such rotations.
%  rotangles(2,sphere_rot), used below, rotates all inclusions to be 
%  parallel to the x-axis instead of the z-axis (z -> x, x -> -z, y -> y) or, 
%              equivalently, Euler angles= (0,pi/2,0).
%

    %
    % set geometry
    %
    NumSpheres = 2; % two pseudoparticles
    center     = zeros(3, NumSpheres);
    x          =-0.5*SepDist(1);
    y          =-0.5*SepDist(2);
    z          =-0.5*SepDist(3);
    for i = 1:NumSpheres
        center(1:3, i) = [x,y,z];
        x = x + SepDist(1);
        y = y + SepDist(2);
        z = z + SepDist(3);
    end;

    %
    % set radius and type info
    %
    radius      = zeros(1,NumSpheres);
    radius(:)   = SphereRadius ; 
    sphere_type = zeros(1,NumSpheres);
    for i = 1:NumSpheres
        sphere_type(1,i) = SphereType;
    end

    %
    % If type equals 2, set eps and mu for each dielectric sphere.
    %
    sphere_eps = zeros(1,NumSpheres)+1i*zeros(1,NumSpheres);
    sphere_cmu = zeros(1,NumSpheres)+1i*zeros(1,NumSpheres);

    for i = 1:NumSpheres
%         sphere_eps(1,i) = 1;
%         sphere_cmu(1,i) = 1;
        sphere_eps(1,i) = SphereEpr;
        sphere_cmu(1,i) = SphereMur;
    end

    %
    %  Change orientation of inclusions, if desired.
    %
    sphere_rot = zeros(3,NumSpheres);  
    for i = 1:NumSpheres
        sphere_rot(1:3,i) = [0,0,0];      
    end
    sphere_rot = rotangles(1,sphere_rot); % 1 = do nothing
    
    
    %
    % save results to a single variable
    %
    A.nspheres = NumSpheres;
    A.center   = center;
    A.radius   = radius;
    A.type     = sphere_type;
    A.eps      = sphere_eps;
    A.cmu      = sphere_cmu;
    A.rot      = sphere_rot;
    A.ncoefs   = (nterms+1)*(2*nterms+1);
    
    A.eps0     = 1;
    A.cmu0     = 1;
    
    % output results
    varargout{1} = A;
    
end