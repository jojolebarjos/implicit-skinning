
% Load some spheres
spheres = cell(6, 1);
spheres{1}.center = [3, 2, 1];
spheres{1}.radius = 1;
spheres{2}.center = [-2, 0, 1.5];
spheres{2}.radius = 2;
spheres{3}.center = [1, 1, 0];
spheres{3}.radius = 0.5;
spheres{4}.center = [0, 0, 3];
spheres{4}.radius = 0.8;
spheres{5}.center = [0, 0, -2];
spheres{5}.radius = 1.5;
spheres{6}.center = [1, -3, 1];
spheres{6}.radius = 1.2;
blocks = cell(3, 1);
blocks{1} = [1];
blocks{2} = [2, 3];
blocks{3} = [4, 5, 6];

% Show SDF
[X, Y, Z] = meshgrid(-5 : 0.25 : 5, -5 : 0.25 : 5, -5 : 0.25 : 5);
V = reshape(sphere_apply(spheres, blocks, [X(:), Y(:), Z(:)]), size(X));
p = patch(isosurface(X, Y, Z, V, 0));
isonormals(X, Y, Z, V, p);
p.FaceColor = 'red';
p.FaceAlpha = 0.5;
p.EdgeColor = 'none';
view([-90, 0]);
camlight;
view([90, 0]);
camlight;
axis equal;
grid off;
lighting gouraud;
axis off;

