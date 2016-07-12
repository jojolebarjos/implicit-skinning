% Register gptoolbox
addpath(genpath('external'));

% Load mesh
load('mesh/neutral.mat');

% Define some transform
reference = {};
[reference.vertices, reference.faces] = load_mesh('mesh/3.obj');
load('mesh/centers3.mat');
centers = zeros(length(centers3), 3);
for i = 1 : length(centers3)
    centers(i, :) = [centers3{i}(2), centers3{i}(3), centers3{i}(1)];
end
transforms = bone_transforms(mesh, centers);

% Use linear skinning
linear = skin_linear(mesh, transforms);
trimesh(reference.faces, reference.vertices(:, 1), reference.vertices(:, 2), reference.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.0, 0.8, 0.0], 'FaceAlpha', 0.2);
hold on;
trimesh(linear.faces, linear.vertices(:, 1), linear.vertices(:, 2), linear.vertices(:, 3), 'EdgeColor', [0.4, 0.25, 0.0], 'EdgeAlpha', 0.2, 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 1.0);
% quiver3(linear.vertices(:, 1), linear.vertices(:, 2), linear.vertices(:, 3), linear.normals(:, 1), linear.normals(:, 2), linear.normals(:, 3), 'Color', [0.8, 0.5, 0.0]);
hold off;
view([-90, 0]);
camlight;
view([90, 0]);
camlight;
axis equal;
grid off;
lighting gouraud;
axis off;
title('Linear skinning');
pause;

% Use dual quaternion skinning
dualquat = skin_dualquat(mesh, transforms);
v = get(gca, 'view');
trimesh(reference.faces, reference.vertices(:, 1), reference.vertices(:, 2), reference.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.0, 0.8, 0.0], 'FaceAlpha', 0.2);
hold on;
trimesh(dualquat.faces, dualquat.vertices(:, 1), dualquat.vertices(:, 2), dualquat.vertices(:, 3), 'EdgeColor', [0.4, 0.25, 0.0], 'EdgeAlpha', 0.2, 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 1.0);
% quiver3(dualquat.vertices(:, 1), dualquat.vertices(:, 2), dualquat.vertices(:, 3), dualquat.normals(:, 1), dualquat.normals(:, 2), dualquat.normals(:, 3), 'Color', [0.8, 0.5, 0.0]);
hold off;
view([-90, 0]);
camlight;
view([90, 0]);
camlight;
axis equal;
grid off;
lighting gouraud;
axis off;
title('Dual quaternion skinning');
set(gca, 'view', v);
pause;

% Use implicit skinning
implicit = skin_implicit(mesh, transforms, 1, 'hrbf');
% implicit = skin_implicit(mesh, transforms, 1, 'spheres');
v = get(gca, 'view');
trimesh(reference.faces, reference.vertices(:, 1), reference.vertices(:, 2), reference.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.0, 0.8, 0.0], 'FaceAlpha', 0.2);
hold on;
trimesh(implicit.faces, implicit.vertices(:, 1), implicit.vertices(:, 2), implicit.vertices(:, 3), 'EdgeColor', [0.4, 0.25, 0.0], 'EdgeAlpha', 0.2, 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 1.0);
% quiver3(implicit.vertices(:, 1), implicit.vertices(:, 2), implicit.vertices(:, 3), implicit.normals(:, 1), implicit.normals(:, 2), implicit.normals(:, 3), 'Color', [0.8, 0.5, 0.0]);
hold off;
view([-90, 0]);
camlight;
view([90, 0]);
camlight;
axis equal;
grid off;
lighting gouraud;
axis off;
title('Implicit skinning');
set(gca, 'view', v);
