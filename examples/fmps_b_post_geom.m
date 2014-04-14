
id=1;

% big sphere
plot_spheres(id,'sphere320.a.tri',1,[0 0 0]',radius0);

plot_spheres(id,'sphere180.a.tri',nspheres,center,radius);

plot_inclusions(id,'2ellipsoids-25x25x75-sep40-draft.a.tri',nspheres,center,radius,sphere_rot);
