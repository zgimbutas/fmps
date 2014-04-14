%
% script runs fmps software
%
% modified by Maxim Ignatenko (KSU, 06/2011)
%
function varargout = FMPS_run(A, wavelength, PropDir, Polarization)

      % read input data
      nnodes     = A.nnodes;
      nphi       = A.nphi;
      ntheta     = A.ntheta;
      rnodes     = A.rnodes;
      weights    = A.weights;
      nterms     = A.nterms; % order of decomposition
      center     = A.center; % centers of pseudoparticles
      radius     = A.radius; % radii of pseudoparticles
      eps0       = A.eps0;
      cmu0       = A.cmu0;
      nspheres   = length(radius);
      
      omega      = 2*pi/wavelength ;      % angular frequaency
      zk         = omega*sqrt(eps0*cmu0); % wavenumber

      fprintf(['wavelength = ',num2str(wavelength)]) ; 
      
      % Set up remaining solver parameters
      A.omega    = omega;
      A.zk       = zk;

    %  STEP 4A: Define the incoming field.
    %  The incoming plane wave is specified by the direction vector (kvec) and
    %  the corresponding E polarization vector (epol).
      kvec       = zk*PropDir;
      epol       = Polarization;

    %  STEP 4B: get incoming E and H fields on all FMPS sphere boundaries

      evecs  = zeros(3,nnodes,nspheres) + 1i*zeros(3,nnodes,nspheres);
      hvecs  = zeros(3,nnodes,nspheres) + 1i*zeros(3,nnodes,nspheres);

      for j = 1:nspheres
        ntargets = nphi*ntheta;
        targets  = repmat(center(:,j), 1, nnodes) + rnodes*radius(j);

        [evecs(:,:,j),hvecs(:,:,j)] = ...
          em3d_planearb_targeval(kvec, epol, ntargets, targets);
      end

    %  STEP 4C: convert incoming E and H fields to local vector 
    %           spherical harmonic expansions.
      ncoefs  = ( nterms + 1 ) * ( 2*nterms + 1 );
      aimpole = zeros(ncoefs, nspheres);
      bimpole = zeros(ncoefs, nspheres);

      for j = 1:nspheres
        [aimpole(:,j),bimpole(:,j)]=...
          em3ehformta(zk, center(:,j), radius(j), ...
          rnodes, weights, nphi, ntheta,  ...
          evecs(:,:,j), hvecs(:,:,j), nterms);
      end

    %
    %  Matlab + Fortran 90 solver: Extract relevant scattering matrix from a file in the
    %  directory
    %
    %  ../data/scat.2ellipsoids-25x25x75-sep40-gold-draft
    %

      [A.raa,A.rab,A.rba,A.rbb,A.nterms_r,ier]=rmatr_file('fort.19');

    %
    %  Construct the right hand side for scattering problem.
    %  by applying the scattering/reflection matrix to incoming data.
    %
      [arhs,brhs]=em3d_multa_r(center,radius,aimpole,bimpole,nterms,A);

    %
    %  Call the solver em3d_multa_fmps with problem now set up.
    %
    %  fprintf('FMPS solver for the Maxwell equation in R^3, Matlab + Fortran 90')
    %  tic

      rhs = reshape([arhs brhs], ncoefs * nspheres*2, 1);
      sol = gmres_simple(@(x) em3d_multa_fmps(A,x), rhs, 1e-3, 20);

      sol0 = reshape(sol, ncoefs, nspheres,2);
      aompole = sol0(:,:,1);
      bompole = sol0(:,:,2);

    %  time_gmres=toc

    %
    %  Generate incoming multipole expansions from all (now known) scattering expansions.
    %
      [asmpole,bsmpole] = em3d_multa_mptaf90(nspheres,nterms,ncoefs,...
          omega,eps0,cmu0,center,radius,...
          aompole,bompole,rnodes,weights,nphi,ntheta);


      varargout{1} = aimpole;
      varargout{2} = bimpole;
      varargout{3} = aompole;
      varargout{4} = bompole;
      varargout{5} = asmpole;
      varargout{6} = bsmpole;
      varargout{7} = A;
end
