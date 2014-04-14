%
% script calculates x-sections based on the results of fmps calculations
%
% modified by Maxim Ignatenko (KSU, 06/2011)
%
function varargout = FMPS_post_Q(A,aimpole,bimpole,aompole,bompole,asmpole,bsmpole)

%     eps0     = 8.8541878e-12 ;
%     cmu0     = 4*pi*1e-7 ; 
%     clight   = 2.99792458e8 ;
    
    nnodes   = A.nnodes;
    rnodes   = A.rnodes;
    nphi     = A.nphi;
    ntheta   = A.ntheta;
    weights  = A.weights;
    nterms   = A.nterms;       % order of decompositions
    zk       = A.zk;           % wavenumber of incident wave
%     wavelength = A.omega/2/pi; % incident wavelength [nm]
    nspheres = A.nspheres;     % number of pseusoparticles
    center   = A.center;       % positions of the centers of pseudoparticles
    radius   = A.radius;       % radii of the pseudoparticles

%
%  POSTPROCESSING
%
%  Evaluate E and H fields on all spheres
%
    eivecs = zeros(3,nnodes,nspheres) + 1i*zeros(3,nnodes,nspheres);
    hivecs = zeros(3,nnodes,nspheres) + 1i*zeros(3,nnodes,nspheres);
    esvecs = zeros(3,nnodes,nspheres) + 1i*zeros(3,nnodes,nspheres);
    hsvecs = zeros(3,nnodes,nspheres) + 1i*zeros(3,nnodes,nspheres);
    eovecs = zeros(3,nnodes,nspheres) + 1i*zeros(3,nnodes,nspheres);
    hovecs = zeros(3,nnodes,nspheres) + 1i*zeros(3,nnodes,nspheres);

    for i=1:nspheres
        [evecs,hvecs]=...
          em3taevaleh(zk,center(:,i),aimpole(:,i),bimpole(:,i),nterms, ...
          center(:,i),radius(i),rnodes,weights,nphi,ntheta);
        eivecs(:,:,i) = evecs; % E incident
        hivecs(:,:,i) = hvecs; % H incident

        [evecs,hvecs]=...
          em3taevaleh(zk,center(:,i),asmpole(:,i),bsmpole(:,i),nterms, ...
          center(:,i),radius(i),rnodes,weights,nphi,ntheta);
        esvecs(:,:,i) = evecs;
        hsvecs(:,:,i) = hvecs;

        [evecs,hvecs]=...
          em3mpevaleh(zk,A.center(:,i),aompole(:,i),bompole(:,i),nterms, ...
          center(:,i),radius(i),rnodes,weights,nphi,ntheta);
        eovecs(:,:,i) = evecs;
        hovecs(:,:,i) = hvecs;
    end;

% total field 
    etvecs = eivecs + esvecs + eovecs;
    htvecs = hivecs + hsvecs + hovecs;

% scattered field
    eavecs = esvecs + eovecs;
    havecs = hsvecs + hovecs;

%
%  Evaluate absorption, scattering, and extinction cross sections
%
 
    qabs = zeros(1,nspheres) + 1i*zeros(1,nspheres);
    qsca = zeros(1,nspheres) + 1i*zeros(1,nspheres);
    qext = zeros(1,nspheres) + 1i*zeros(1,nspheres);

    for i=1:nspheres
        ptvecs  = cross(etvecs(:,:,i),conj(htvecs(:,:,i)));
        qabs(i) = weights  * sum(ptvecs .* rnodes,1)' / pi;

        pavecs = cross(eavecs(:,:,i),conj(havecs(:,:,i)));
        qsca(i) = weights  * sum(pavecs .* rnodes,1)' / pi;

        pxvecs = cross(eavecs(:,:,i),conj(hivecs(:,:,i))) ...
          + cross(eivecs(:,:,i),conj(havecs(:,:,i)));
        qext(i) = -weights  * sum(pxvecs .* rnodes,1)' / pi;
    end;

%     qabs_total(ifreq) = sum(qabs * pi .* radius.^2);
%     qsca_total(ifreq) = sum(qsca * pi .* radius.^2);
%     qext_total(ifreq) = sum(qext * pi .* radius.^2);
    qabs_total = sum(qabs * pi .* radius.^2);
    qsca_total = sum(qsca * pi .* radius.^2);
    qext_total = sum(qext * pi .* radius.^2);

    qabs_total = qabs_total * (1e-18/2/376.73);
    qsca_total = qsca_total * (1e-18/2/376.73);
    qext_total = qext_total * (1e-18/2/376.73);   

    varargout{1} = qabs_total; % cross-sections
    varargout{2} = qsca_total;
    varargout{3} = qext_total;

end