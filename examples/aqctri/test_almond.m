
nsource = 1000;
source = zeros(3,nsource);

  theta=rand(1,nsource)*pi;
  phi=rand(1,nsource)*2*pi;
  source(1,:)=.5*cos(phi).*sin(theta);
  source(2,:)=.5*sin(phi).*sin(theta);
  source(3,:)=.5*cos(theta);

points=sphere2almond(source);

scatter3(points(1,:),points(2,:),points(3,:))
