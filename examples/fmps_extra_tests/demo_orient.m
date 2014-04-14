% test rotangles itype=2 orientation

% x-axis orientation (z -> x, x -> -z, y -> y)

rotmat=em3orient([0 pi/2 0]',1);

rotmat*[1 0 0]'
rotmat*[0 1 0]'
rotmat*[0 0 1]'
