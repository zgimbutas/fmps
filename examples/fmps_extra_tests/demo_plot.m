
[nspheres,sphere_xyz,sphere_r]=demo_geometry(7);

%
%
%  Change orientation of inclusions, if desired.
%
  
sphere_rot = zeros(3,nspheres);
for i=1:nspheres
    sphere_rot(1:3,i)=[0,0,0];
end

sphere_rot=rotangles(2,sphere_rot);

fprintf('nspheres = %d\n',nspheres)


id=1;

% big sphere
%plot_spheres(id,'sphere320.a.tri',1,[0 0 0]',2480/4);

plot_spheres(id,'sphere180.a.tri',nspheres,sphere_xyz,sphere_r);

plot_inclusions(id,'2ellipsoids-25x25x75-sep40-draft.a.tri',nspheres,sphere_xyz,sphere_r,sphere_rot);
