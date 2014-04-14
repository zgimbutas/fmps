function points_m = sphere2dogive(points)
% Double ogive
%
%  Translated from Felipe Vico's Fortran tests.
%
%  This program allows to generate the geometries defined in [1] 
%  and commonly used to test commercial and scientific EM codes 
%  in the IEEE community:
%
%  http://www.efieldsolutions.com/example_rcs_almond.pdf
%  http://www.feko.info/applications/white-papers/
%        rcs-benchmarking-of-generic-shapes
%  http://www.cst.com/Content/Applications/Article/article.aspx?id=133
%
%  [1] A.C. Woo, H.T. Wang, M.J. Schuh, M.L. Sanders,
%  'Benchmark radar targets for the validation of computational 
%  electromagnetic programs'
%  IEEE Antennas and Propagation Magazine, Vol.35 No. 1 February 1993

  m=size(points,2);
  points_m=zeros(3,m);

  for i=1:m
    x=points(1,i)*6;
    y=points(2,i);
    z=points(3,i);
    
    if (x > 0),
      c4=13;
      c5=13;
      c6=13;
      alpha=0.9230760;
    else
      c4=3.1957816;
      c5=3.4408025;
      c6=3.1957816;
      alpha=0.6870875;
    end

    K=(y/c4)^2+(z/c6)^2;
    MM=4*K*alpha^2-4*(alpha^2-1)*(K+(x/c5)^2);
    p2=(-2*sqrt(K)*alpha+sqrt(MM))/(2*(K+(x/c5)^2));
    x=x*p2;
    y=y*p2;
    z=z*p2;
       
    points_m(1,i)=x;
    points_m(2,i)=y;
    points_m(3,i)=z;
  end

