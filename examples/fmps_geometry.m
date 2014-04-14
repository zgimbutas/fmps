function [nspheres,sphere_xyz,sphere_r]=fmps_geometry(itype)
%
% Generate a distribution of spheres.
%
% itype = 1..17
%
%    1: one sphere, r=50
%    2: two spheres, r=50, sep=200
%    3: three spheres, r=50, sep=200
%    4: four spheres, 2x2, r=50, sep=200
%    5: hexagonal grid, 9x9x2, r=50, sep=105
%    6: regular grid, 5x5x2, r=50, sep=140
%    7: regular grid, 15x15, r=50, sep=140
%    8: regular grid + perturbations, 23x23, r=50, sep=140
%    9: regular grid, 21x21x2, r=50, sep=140
%   10: regular grid + perturbations, 21x21x2, r=50, sep=120
%   11: regular grid, 21x21x2, r=50, sep=120
%   12: regular grid + perturbations, 21x21x2, r=65, sep=156
%   13: regular grid, 21x21x4, r=50, see code
%   14: regular grid, 21x21x4, r=50, sep_xy=240, sep_z=145
%   15: regular grid, 21x21x4, r=65, sep=156, randomized (very slow)
%   16: regular grid, 21x21x4, r=50, sep_xy=120, sep_z=120
%   17: regular grid, 21x21x4, r=50, sep=120, randomized (very slow)
%   18: regular grid, 5x5x5, r=50, sep=140
%   19: regular grid, 9x9x9, r=50, sep=140
%   20: regular grid, 11x11x11, r=50, sep=140
%   21: hexagonal grid, 23x23x1, r=50, sep=105
%   22: hexagonal grid, 41x41x1, r=50, sep=105
%   23: regular grid, 5x5x1, r=50, sep=140
%   24: regular grid, 11x11x1, r=50, sep=140
%   25: hexagonal grid, 9x9x1, r=50, sep=105
%

if( itype == 1 )
%
%  one particle
%

  nspheres=1;
  
  sphere_xyz = zeros(3,nspheres);
  sphere_xyz(1:3,1) = [0,0,0];
  
  sphere_r = zeros(1,nspheres);
  sphere_r(1) = 50;
  

elseif( itype == 2 )
%
%  two particles
%

  nspheres=2;
  
  sphere_xyz = zeros(3,nspheres);
  sphere_xyz(1:3,1) = [-100,0,0];
  sphere_xyz(1:3,2) = [ 100,0,0];
  
  sphere_r = zeros(1,nspheres);
  sphere_r(1) = 50;
  sphere_r(2) = 50;
  

elseif( itype == 3 )
%
%  three particles
%

  nspheres=3;
  
  sphere_xyz = zeros(3,nspheres);
  sphere_xyz(1:3,1) = [0,0,0];
  sphere_xyz(1:3,2) = [0,200,0];
  sphere_xyz(1:3,3) = [200,0,0];
  
  sphere_r = zeros(1,nspheres);
  sphere_r(1) = 50;
  sphere_r(2) = 50;
  sphere_r(3) = 50;
  

elseif( itype == 4 )
%
%  four particles
%

  nspheres=4;
  
  sphere_xyz = zeros(3,nspheres);
  sphere_xyz(1:3,1) = [0,-100,0];
  sphere_xyz(1:3,2) = [0,+100,0];
  sphere_xyz(1:3,3) = [-100,0,0];
  sphere_xyz(1:3,4) = [+100,0,0];
  
  sphere_r = zeros(1,nspheres);
  sphere_r(1) = 50;
  sphere_r(2) = 50;
  sphere_r(3) = 50;
  sphere_r(4) = 50;
  

elseif( itype == 5 )
%
%  hexagonal regular grid 
%

  radius=50;
  rsep=radius*1.05;
  
  ngridx=5;
  ngridy=5;
  ngridz=2;

  nspheres=(2*ngridx-1)*(2*ngridy-1)*(ngridz);

  sphere_xyz = zeros(3,nspheres);
  sphere_r = zeros(1,nspheres);
  

  kk=0;

  for ii=-ngridx+1:ngridx-1
  for jj=-ngridy+1:ngridy-1
  for ll=0:ngridz-1
    kk=kk+1;
    if( mod(jj,2) == 0 ) 
        sphere_xyz(1,kk)=+2*rsep*(ii);
    else
        sphere_xyz(1,kk)=+2*rsep*(ii+.5);
    end
    if( mod(ll,2) == 0 ) 
      sphere_xyz(2,kk)=+2*rsep*(jj)*sqrt(3.0)/2;
    else
      sphere_xyz(2,kk)=+2*rsep*(jj+1)*sqrt(3.0)/2;
    end
    sphere_xyz(3,kk)=+2*rsep*(ll)*sqrt(3.0)/2;
    sphere_r(kk)=radius;
  end
  end
  end


elseif( itype == 6 )
%
%  regular grid 
%

  radius=50;
  rsep=radius*1.4;
  
  ngridx=3;
  ngridy=3;
  ngridz=2;

  nspheres=(2*ngridx-1)*(2*ngridy-1)*(ngridz);

  sphere_xyz = zeros(3,nspheres);
  sphere_r = zeros(1,nspheres);
  
  kk=0;

  for ii=-ngridx+1:ngridx-1
  for jj=-ngridy+1:ngridy-1
  for ll=0:ngridz-1
    kk=kk+1;
    sphere_xyz(1,kk)=+2*rsep*(ii);
    sphere_xyz(2,kk)=+2*rsep*(jj);
    sphere_xyz(3,kk)=+2*rsep*(ll);
    sphere_r(kk)=radius;
  end
  end
  end


elseif( itype == 7 )
%
%  regular grid 
%

  radius=50;
  rsep=radius*1.4;
  
  ngridx=8;
  ngridy=8;
  ngridz=1;

  nspheres=(2*ngridx-1)*(2*ngridy-1)*(ngridz);

  sphere_xyz = zeros(3,nspheres);
  sphere_r = zeros(1,nspheres);
  
  kk=0;

  for ii=-ngridx+1:ngridx-1
  for jj=-ngridy+1:ngridy-1
  for ll=0:ngridz-1
    kk=kk+1;
    sphere_xyz(1,kk)=+2*rsep*(ii);
    sphere_xyz(2,kk)=+2*rsep*(jj);
    sphere_xyz(3,kk)=+2*rsep*(ll);
    sphere_r(kk)=radius;
  end
  end
  end


elseif( itype == 8 )
%
%  regular grid + random perturbations 
%

  radius=50;
  rsep=radius*1.4;
  
  ngridx=12;
  ngridy=12;
  ngridz=1;

  nspheres=(2*ngridx-1)*(2*ngridy-1)*(ngridz);

  sphere_xyz = zeros(3,nspheres);
  sphere_r = zeros(1,nspheres);
  
  kk=0;

  for ii=-ngridx+1:ngridx-1
  for jj=-ngridy+1:ngridy-1
  for ll=0:ngridz-1
    kk=kk+1;
    sphere_xyz(1,kk)=+2*rsep*(ii) + 2*20*(rand-0.5);
    sphere_xyz(2,kk)=+2*rsep*(jj) + 2*20*(rand-0.5);
    sphere_xyz(3,kk)=+2*rsep*(ll);
    sphere_r(kk)=radius;
  end
  end
  end


elseif( itype == 9 )
%
%  regular grid 
%

  radius=50;
  rsep=radius*1.4;
  
  ngridx=11;
  ngridy=11;
  ngridz=2;

  nspheres=(2*ngridx-1)*(2*ngridy-1)*(ngridz);

  sphere_xyz = zeros(3,nspheres);
  sphere_r = zeros(1,nspheres);
  
  kk=0;

  for ii=-ngridx+1:ngridx-1
  for jj=-ngridy+1:ngridy-1
  for ll=0:ngridz-1
    kk=kk+1;
    sphere_xyz(1,kk)=+2*rsep*(ii);
    sphere_xyz(2,kk)=+2*rsep*(jj);
    sphere_xyz(3,kk)=+2*rsep*(ll);
    sphere_r(kk)=radius;
  end
  end
  end


elseif( itype == 10 )
%
%  regular grid + random perturbations 
%

  radius=50;
  rsep=radius*1.2;
  
  ngridx=11;
  ngridy=11;
  ngridz=2;

  nspheres=(2*ngridx-1)*(2*ngridy-1)*(ngridz);

  sphere_xyz = zeros(3,nspheres);
  sphere_r = zeros(1,nspheres);
  
  kk=0;

  for ii=-ngridx+1:ngridx-1
  for jj=-ngridy+1:ngridy-1
  for ll=0:ngridz-1
    kk=kk+1;
    sphere_xyz(1,kk)=+2*rsep*(ii) + 2*10*(rand-0.5);
    sphere_xyz(2,kk)=+2*rsep*(jj) + 2*10*(rand-0.5);
    sphere_xyz(3,kk)=+2*rsep*(ll);
    sphere_r(kk)=radius;
  end
  end
  end


elseif( itype == 11 )
%
%  regular grid 
%

  radius=50;
  rsep=radius*1.2;
  
  ngridx=11;
  ngridy=11;
  ngridz=2;

  nspheres=(2*ngridx-1)*(2*ngridy-1)*(ngridz);

  sphere_xyz = zeros(3,nspheres);
  sphere_r = zeros(1,nspheres);
  
  kk=0;

  for ii=-ngridx+1:ngridx-1
  for jj=-ngridy+1:ngridy-1
  for ll=0:ngridz-1
    kk=kk+1;
    sphere_xyz(1,kk)=+2*rsep*(ii);
    sphere_xyz(2,kk)=+2*rsep*(jj);
    sphere_xyz(3,kk)=+2*rsep*(ll);
    sphere_r(kk)=radius;
  end
  end
  end


elseif( itype == 12 )
%
%  regular grid + random perturbations 
%

  radius=65;
  rsep=radius*1.2;
  
  ngridx=11;
  ngridy=11;
  ngridz=2;

  nspheres=(2*ngridx-1)*(2*ngridy-1)*(ngridz);

  sphere_xyz = zeros(3,nspheres);
  sphere_r = zeros(1,nspheres);
  
  kk=0;

  for ii=-ngridx+1:ngridx-1
  for jj=-ngridy+1:ngridy-1
  for ll=0:ngridz-1
    kk=kk+1;
    sphere_xyz(1,kk)=+2*rsep*(ii) + 2*20*(rand-0.5);
    sphere_xyz(2,kk)=+2*rsep*(jj) + 2*20*(rand-0.5);
    sphere_xyz(3,kk)=+2*rsep*(ll);
    sphere_r(kk)=radius;
  end
  end
  end


elseif( itype == 13 )
%
%  regular grid 
%

  radius=50;
  rsep=radius*1.2;
  
  ngridx=11;
  ngridy=11;
  ngridz=4;

  nspheres=(2*ngridx-1)*(2*ngridy-1)*(ngridz);

  sphere_xyz = zeros(3,nspheres);
  sphere_r = zeros(1,nspheres);
  
  kk=0;

  for ii=-ngridx+1:ngridx-1
  for jj=-ngridy+1:ngridy-1
  for ll=0:ngridz-1
    kk=kk+1;
    sphere_xyz(1,kk)=2*(100*1.2d0)*(ii);
    sphere_xyz(2,kk)=2*(100*1.2d0)*(jj);
    if( ll == 0  )  sphere_xyz(3,kk)=-72.5; end
    if( ll == 1  )  sphere_xyz(3,kk)=+72.5; end
    if( ll == 2  )  sphere_xyz(3,kk)=-72.5 + 2*(100*1.2); end
    if( ll == 3  )  sphere_xyz(3,kk)=+72.5 + 2*(100*1.2); end
    sphere_r(kk)=radius;
  end
  end
  end


elseif( itype == 14 )
%
%  regular grid 
%

  radius=50;
  rsep=72.5;
  
  ngridx=11;
  ngridy=11;
  ngridz=4;

  nspheres=(2*ngridx-1)*(2*ngridy-1)*(ngridz);

  sphere_xyz = zeros(3,nspheres);
  sphere_r = zeros(1,nspheres);
  
  kk=0;

  for ii=-ngridx+1:ngridx-1
  for jj=-ngridy+1:ngridy-1
  for ll=0:ngridz-1
    kk=kk+1;
    sphere_xyz(1,kk)=2*(100*1.2d0)*(ii);
    sphere_xyz(2,kk)=2*(100*1.2d0)*(jj);
    if( ll == 0  )  sphere_xyz(3,kk)=-72.5; end
    if( ll == 1  )  sphere_xyz(3,kk)=+72.5*1; end
    if( ll == 2  )  sphere_xyz(3,kk)=+72.5*3; end
    if( ll == 3  )  sphere_xyz(3,kk)=+72.5*5; end
    sphere_r(kk)=radius;
  end
  end
  end


elseif( itype == 15 )
%
%  regular grid + randomized locations (very slow under matlab)
%

  radius=65;
  rsep=radius*1.2;
  
  ngridx=11;
  ngridy=11;
  ngridz=2;

  nspheres=(2*ngridx-1)*(2*ngridy-1)*(ngridz);

  sphere_xyz = zeros(3,nspheres);
  sphere_r = zeros(1,nspheres);
  
  kk=0;

  for ii=-ngridx+1:ngridx-1
  for jj=-ngridy+1:ngridy-1
  for ll=0:ngridz-1
    kk=kk+1;
    sphere_xyz(1,kk)=+2*rsep*(ii) + 2*10*(rand-0.5);
    sphere_xyz(2,kk)=+2*rsep*(jj) + 2*10*(rand-0.5);
    sphere_xyz(3,kk)=+2*rsep*(ll);
    sphere_r(kk)=radius;
  end
  end
  end

  center = zeros(3,1)
  for niter=1:100
  for jj=1:nspheres

    center(1) = sphere_xyz(1,jj) + 10*2*(rand-0.5);
    center(2) = sphere_xyz(2,jj) + 10*2*(rand-0.5);
    center(3) = sphere_xyz(3,jj); 
        
    ifinter=0;

    for ii=1:nspheres       
      if( ii == jj ) continue; end
      d=(center(1)-sphere_xyz(1,ii))^2;
      d=d+(center(2)-sphere_xyz(2,ii))^2;
      d=d+(center(3)-sphere_xyz(3,ii))^2;
      d=sqrt(d);
      if( d < (sphere_r(jj)+sphere_r(ii)) ) 
      ifinter=1;
      continue
      end
    end

    if( ifinter == 0 ) 
      sphere_xyz(1,jj)=center(1);
      sphere_xyz(2,jj)=center(2);
      sphere_xyz(3,jj)=center(3);
    end
        
  end
  end


elseif( itype == 16 )
%
%  regular grid 
%

  radius=50;
  rsep=60;
  
  ngridx=11;
  ngridy=11;
  ngridz=4;

  nspheres=(2*ngridx-1)*(2*ngridy-1)*(ngridz);

  sphere_xyz = zeros(3,nspheres);
  sphere_r = zeros(1,nspheres);
  
  kk=0;

  for ii=-ngridx+1:ngridx-1
  for jj=-ngridy+1:ngridy-1
  for ll=0:ngridz-1
    kk=kk+1;
    sphere_xyz(1,kk)=2*rsep*(ii);
    sphere_xyz(2,kk)=2*rsep*(jj);
    if( ll == 0  )  sphere_xyz(3,kk)=-60; end
    if( ll == 1  )  sphere_xyz(3,kk)=+60*1; end
    if( ll == 2  )  sphere_xyz(3,kk)=+60*3; end
    if( ll == 3  )  sphere_xyz(3,kk)=+60*5; end
    sphere_r(kk)=radius;
  end
  end
  end


elseif( itype == 17 )
%
%  regular grid + randomized locations (very slow under matlab)
%  dense packing
%

  radius=50;
  rsep=radius*1.05;
  
  ngridx=11;
  ngridy=11;
  ngridz=2;

  nspheres=(2*ngridx-1)*(2*ngridy-1)*(ngridz);

  sphere_xyz = zeros(3,nspheres);
  sphere_r = zeros(1,nspheres);
  
  kk=0;

  for ii=-ngridx+1:ngridx-1
  for jj=-ngridy+1:ngridy-1
  for ll=0:ngridz-1
    kk=kk+1;
    sphere_xyz(1,kk)=+2*rsep*(ii) + 2*2.5*(rand-0.5);
    sphere_xyz(2,kk)=+2*rsep*(jj) + 2*2.5*(rand-0.5);
    sphere_xyz(3,kk)=+2*rsep*(ll);
    sphere_r(kk)=radius;
  end
  end
  end

  center = zeros(3,1)
  for niter=1:100
  for jj=1:nspheres

    center(1) = sphere_xyz(1,jj) + 2*2.5*(rand-0.5);
    center(2) = sphere_xyz(2,jj) + 2*2.5*(rand-0.5);
    center(3) = sphere_xyz(3,jj); 
        
    ifinter=0;

    for ii=1:nspheres       
      if( ii == jj ) 
          continue; 
      end
      d=(center(1)-sphere_xyz(1,ii))^2;
      d=d+(center(2)-sphere_xyz(2,ii))^2;
      d=d+(center(3)-sphere_xyz(3,ii))^2;
      d=sqrt(d);
      if( d < (sphere_r(jj)+sphere_r(ii)) ) 
      ifinter=1;
      continue
      end
    end

    if( ifinter == 0 ) 
      sphere_xyz(1,jj)=center(1);
      sphere_xyz(2,jj)=center(2);
      sphere_xyz(3,jj)=center(3);
    end
        
  end
  end


elseif( itype == 18 )
%
%  regular grid 
%

  radius=50;
  rsep=radius*1.4;
  
  ngridx=3;
  ngridy=3;
  ngridz=3;

  nspheres=(2*ngridx-1)*(2*ngridy-1)*(2*ngridz-1);

  sphere_xyz = zeros(3,nspheres);
  sphere_r = zeros(1,nspheres);
  
  kk=0;

  for ii=-ngridx+1:ngridx-1
  for jj=-ngridy+1:ngridy-1
  for ll=-ngridz+1:ngridz-1
    kk=kk+1;
    sphere_xyz(1,kk)=+2*rsep*(ii);
    sphere_xyz(2,kk)=+2*rsep*(jj);
    sphere_xyz(3,kk)=+2*rsep*(ll);
    sphere_r(kk)=radius;
  end
  end
  end

elseif( itype == 19 )
%
%  regular grid 
%

  radius=50;
  rsep=radius*1.4;
  
  ngridx=5;
  ngridy=5;
  ngridz=5;

  nspheres=(2*ngridx-1)*(2*ngridy-1)*(2*ngridz-1);

  sphere_xyz = zeros(3,nspheres);
  sphere_r = zeros(1,nspheres);
  
  kk=0;

  for ii=-ngridx+1:ngridx-1
  for jj=-ngridy+1:ngridy-1
  for ll=-ngridz+1:ngridz-1
    kk=kk+1;
    sphere_xyz(1,kk)=+2*rsep*(ii);
    sphere_xyz(2,kk)=+2*rsep*(jj);
    sphere_xyz(3,kk)=+2*rsep*(ll);
    sphere_r(kk)=radius;
  end
  end
  end

elseif( itype == 20 )
%
%  regular grid 
%

  radius=50;
  rsep=radius*1.4;
  
  ngridx=6;
  ngridy=6;
  ngridz=6;

  nspheres=(2*ngridx-1)*(2*ngridy-1)*(2*ngridz-1);

  sphere_xyz = zeros(3,nspheres);
  sphere_r = zeros(1,nspheres);
  
  kk=0;

  for ii=-ngridx+1:ngridx-1
  for jj=-ngridy+1:ngridy-1
  for ll=-ngridz+1:ngridz-1
    kk=kk+1;
    sphere_xyz(1,kk)=+2*rsep*(ii);
    sphere_xyz(2,kk)=+2*rsep*(jj);
    sphere_xyz(3,kk)=+2*rsep*(ll);
    sphere_r(kk)=radius;
  end
  end
  end

elseif( itype == 21 )
%
%  hexagonal regular grid 
%

  radius=50;
  rsep=radius*1.05;
  
  ngridx=12;
  ngridy=12;
  ngridz=1;

  nspheres=(2*ngridx-1)*(2*ngridy-1)*(ngridz);

  sphere_xyz = zeros(3,nspheres);
  sphere_r = zeros(1,nspheres);
  

  kk=0;

  for ii=-ngridx+1:ngridx-1
  for jj=-ngridy+1:ngridy-1
  for ll=0:ngridz-1
    kk=kk+1;
    if( mod(jj,2) == 0 ) 
        sphere_xyz(1,kk)=+2*rsep*(ii);
    else
        sphere_xyz(1,kk)=+2*rsep*(ii+.5);
    end
    if( mod(ll,2) == 0 ) 
      sphere_xyz(2,kk)=+2*rsep*(jj)*sqrt(3.0)/2;
    else
      sphere_xyz(2,kk)=+2*rsep*(jj+1)*sqrt(3.0)/2;
    end
    sphere_xyz(3,kk)=+2*rsep*(ll)*sqrt(3.0)/2;
    sphere_r(kk)=radius;
  end
  end
  end


elseif( itype == 22 )
%
%  hexagonal regular grid 
%

  radius=50;
  rsep=radius*1.05;
  
  ngridx=21;
  ngridy=21;
  ngridz=1;

  nspheres=(2*ngridx-1)*(2*ngridy-1)*(ngridz);

  sphere_xyz = zeros(3,nspheres);
  sphere_r = zeros(1,nspheres);
  

  kk=0;

  for ii=-ngridx+1:ngridx-1
  for jj=-ngridy+1:ngridy-1
  for ll=0:ngridz-1
    kk=kk+1;
    if( mod(jj,2) == 0 ) 
        sphere_xyz(1,kk)=+2*rsep*(ii);
    else
        sphere_xyz(1,kk)=+2*rsep*(ii+.5);
    end
    if( mod(ll,2) == 0 ) 
      sphere_xyz(2,kk)=+2*rsep*(jj)*sqrt(3.0)/2;
    else
      sphere_xyz(2,kk)=+2*rsep*(jj+1)*sqrt(3.0)/2;
    end
    sphere_xyz(3,kk)=+2*rsep*(ll)*sqrt(3.0)/2;
    sphere_r(kk)=radius;
  end
  end
  end

elseif( itype == 23 )
%
%  regular grid 
%

  radius=50;
  rsep=radius*1.4;
  
  ngridx=3;
  ngridy=3;
  ngridz=1;

  nspheres=(2*ngridx-1)*(2*ngridy-1)*(2*ngridz-1);

  sphere_xyz = zeros(3,nspheres);
  sphere_r = zeros(1,nspheres);
  
  kk=0;

  for ii=-ngridx+1:ngridx-1
  for jj=-ngridy+1:ngridy-1
  for ll=-ngridz+1:ngridz-1
    kk=kk+1;
    sphere_xyz(1,kk)=+2*rsep*(ii);
    sphere_xyz(2,kk)=+2*rsep*(jj);
    sphere_xyz(3,kk)=+2*rsep*(ll);
    sphere_r(kk)=radius;
  end
  end
  end


elseif( itype == 24 )
%
%  regular grid 
%

  radius=50;
  rsep=radius*1.4;
  
  ngridx=6;
  ngridy=6;
  ngridz=1;

  nspheres=(2*ngridx-1)*(2*ngridy-1)*(2*ngridz-1);

  sphere_xyz = zeros(3,nspheres);
  sphere_r = zeros(1,nspheres);
  
  kk=0;

  for ii=-ngridx+1:ngridx-1
  for jj=-ngridy+1:ngridy-1
  for ll=-ngridz+1:ngridz-1
    kk=kk+1;
    sphere_xyz(1,kk)=+2*rsep*(ii);
    sphere_xyz(2,kk)=+2*rsep*(jj);
    sphere_xyz(3,kk)=+2*rsep*(ll);
    sphere_r(kk)=radius;
  end
  end
  end


elseif( itype == 25 )
%
%  hexagonal regular grid 
%

  radius=50;
  rsep=radius*1.05;
  
  ngridx=5;
  ngridy=5;
  ngridz=1;

  nspheres=(2*ngridx-1)*(2*ngridy-1)*(ngridz);

  sphere_xyz = zeros(3,nspheres);
  sphere_r = zeros(1,nspheres);
  

  kk=0;

  for ii=-ngridx+1:ngridx-1
  for jj=-ngridy+1:ngridy-1
  for ll=0:ngridz-1
    kk=kk+1;
    if( mod(jj,2) == 0 ) 
        sphere_xyz(1,kk)=+2*rsep*(ii);
    else
        sphere_xyz(1,kk)=+2*rsep*(ii+.5);
    end
    if( mod(ll,2) == 0 ) 
      sphere_xyz(2,kk)=+2*rsep*(jj)*sqrt(3.0)/2;
    else
      sphere_xyz(2,kk)=+2*rsep*(jj+1)*sqrt(3.0)/2;
    end
    sphere_xyz(3,kk)=+2*rsep*(ll)*sqrt(3.0)/2;
    sphere_r(kk)=radius;
  end
  end
  end




else

%
%  default, one particle
%

  nspheres=1;
  
  sphere_xyz = zeros(3,nspheres);
  sphere_xyz(1:3,1) = [0,0,0];
  
  sphere_r = zeros(1,nspheres);
  sphere_r(1) = 50;
  
end

