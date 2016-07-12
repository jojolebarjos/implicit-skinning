
% Register gptoolbox
addpath(genpath('external'));

% Load mesh
load('mesh/cylinder.mat');

% Transforms: Bend the cylinder
transforms = cellfun(@double, repmat({eye(4)}, 1, length(mesh.bones)), 'Un', 0);
transforms{2} = matrix_rotation(pi / 1.5, [1, 0, 0], mesh.bones{2}.center);

l1 = [85, 0]; l2 = [95, 0]; vp = [90, 0];

% Use linear skinning
linear = skin_linear(mesh, transforms);
V = linear.vertices; F = linear.faces;
clipping_plane = V(:, 1) > 0;
V(clipping_plane, :) = NaN(nnz(clipping_plane), 3);
trimesh(F, V(:, 1), V(:, 2), V(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 1.0);
view(l1); camlight;
view(l2); camlight;
lighting gouraud;
view(vp);
zoom(2);
axis equal;
grid off;
axis off;
pause;
print('figs/cylinder_linear.png', '-dpng', '-r300');

% Use dual quaternion skinning
dualquat = skin_dualquat(mesh, transforms);
V = dualquat.vertices; F = dualquat.faces;
clipping_plane = V(:, 1) > 0;
V(clipping_plane, :) = NaN(nnz(clipping_plane), 3);
trimesh(F, V(:, 1), V(:, 2), V(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 1.0);
view(l1); camlight;
view(l2); camlight;
lighting gouraud;
view(vp);
zoom(2.8);
axis equal;
grid off;
axis off;
pause;
print('figs/cylinder_dualquat.png', '-dpng', '-r300');

% Use implicit skinning
implicit = skin_implicit(mesh, transforms, 1, 'spheres');
V = implicit.vertices; F = implicit.faces;
clipping_plane = V(:, 1) > 0;
V(clipping_plane, :) = NaN(nnz(clipping_plane), 3);
trimesh(F, V(:, 1), V(:, 2), V(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 1.0);
view(l1); camlight;
view(l2); camlight;
lighting gouraud;
view(vp);
zoom(5);
axis equal;
grid off;
axis off;
pause;
print('figs/cylinder_implicit.png', '-dpng', '-r300');

l1 = [20, -20]; l2 = [95, -20]; vp = [40, -20];
implicit = skin_implicit(mesh, transforms, 1, 'spheres');
V = implicit.vertices; F = implicit.faces;
trimesh(F, V(:, 1), V(:, 2), V(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 1.0, 'SpecularStrength', 0.4);
view(l1); camlight;
view(l2); camlight;
lighting gouraud;
view(vp);
zoom(1.5);
axis equal;
grid off;
axis off;
pause;
print('figs/cylinder_laplacian.png', '-dpng', '-r300');

transforms{2} = matrix_rotation(pi / 1.3, [1, 0, 0], mesh.bones{2}.center);
l1 = [85, 0]; l2 = [0, 0]; vp = [90, 0];

implicit = skin_implicit(mesh, transforms, 1, 'spheres', 'linear');
V = implicit.vertices; F = implicit.faces;
clipping_plane = V(:, 1) > 0;
V(clipping_plane, :) = NaN(nnz(clipping_plane), 3);
trimesh(F, V(:, 1), V(:, 2), V(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 1.0, 'SpecularStrength', 0.6);
view(l1); camlight;
view(l2); camlight;
lighting gouraud;
view(vp);
zoom(1.8);
axis equal;
grid off;
axis off;
pause;
print('figs/cylinder_linear_self.png', '-dpng', '-r300');

implicit = skin_implicit(mesh, transforms, 1, 'spheres', 'dualquat');
V = implicit.vertices; F = implicit.faces;
clipping_plane = V(:, 1) > 0;
V(clipping_plane, :) = NaN(nnz(clipping_plane), 3);
trimesh(F, V(:, 1), V(:, 2), V(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 1.0, 'SpecularStrength', 0.2);
view(l1); camlight;
view([160, 0]); camlight;
view([25, 0]); camlight;
lighting gouraud;
view(vp);
zoom(4.3);
axis equal;
grid off;
axis off;
pause;
print('figs/cylinder_dualquat_self.png', '-dpng', '-r300');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transforms: Twist the cylinder
transforms = cellfun(@double, repmat({eye(4)}, 1, length(mesh.bones)), 'Un', 0);
transforms{2} = matrix_rotation(pi / 1.5, [0, 0, 1], mesh.bones{2}.center);

% Use linear skinning
linear = skin_linear(mesh, transforms);
V = linear.vertices; F = linear.faces;
trimesh(F, V(:, 1), V(:, 2), V(:, 3), 'EdgeColor', 'k', 'EdgeAlpha', 0.2, 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 1.0, 'FaceLighting', 'gouraud', 'EdgeLighting', 'none', 'SpecularStrength', 0.2);
l1 = [75, 0]; view(l1); camlight;
l2 = [105, 0]; view(l2); camlight;
vp = [90, 0]; view(vp);
zoom(2);
axis equal;
grid off;
axis off;
pause;
print('figs/twisted_linear.png', '-dpng', '-r300');

% Use dualquat skinning
dualquat = skin_dualquat(mesh, transforms);
V = dualquat.vertices; F = dualquat.faces;
trimesh(F, V(:, 1), V(:, 2), V(:, 3), 'EdgeColor', 'k', 'EdgeAlpha', 0.2, 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 1.0, 'FaceLighting', 'gouraud', 'EdgeLighting', 'none');
view(l1); camlight;
view(l2); camlight;
view(vp);
zoom(2);
axis equal;
grid off;
axis off;
pause;
print('figs/twisted_dualquat.png', '-dpng', '-r300');

% Use implicit skinning
implicit = skin_implicit(mesh, transforms, 1, 'spheres');
V = implicit.vertices; F = implicit.faces;
trimesh(F, V(:, 1), V(:, 2), V(:, 3), 'EdgeColor', 'k', 'EdgeAlpha', 0.2, 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 1.0, 'FaceLighting', 'gouraud', 'EdgeLighting', 'none');
view(l1); camlight;
view(l2); camlight;
view(vp);
zoom(2);
axis equal;
grid off;
axis off;
pause;
print('figs/twisted_implicit.png', '-dpng', '-r300');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;

% Load mesh
load('mesh/neutral.mat');
mesh = mesh_densify(mesh);
filter = mesh.weights(:, 6) > 0.5 | mesh.weights(:, 7) ~= 0 | mesh.weights(:, 8) ~= 0;

% Define some transform
transforms = cellfun(@double, repmat({eye(4)}, 1, length(mesh.bones)), 'Un', 0);
transforms{7} = matrix_rotation(pi / 3, [0, 0, 1], mesh.bones{7}.center);
transforms{8} = transforms{7} * matrix_rotation(pi / 2, [0, 0, 1], mesh.bones{8}.center);

l1 = [0, 30]; l2 = [180, 30]; vp = [90, 90];

% Use linear skinning
linear = skin_linear(mesh, transforms);
[V, F] = mesh_filter_vertices(linear.vertices, linear.faces, filter);
% clipping_plane = V(:, 3) > 0;
% V(clipping_plane, :) = NaN(nnz(clipping_plane), 3);
trimesh(F, V(:, 1), V(:, 2), V(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 1.0, 'SpecularStrength', 0.4);
view(l1); camlight;
view(l2); camlight;
lighting gouraud;
view(vp);
zoom(1.2);
axis equal;
grid off;
axis off;
pause;
print('figs/index_linear.png', '-dpng', '-r300');

% Use dual quaternion skinning
dualquat = skin_dualquat(mesh, transforms);
[V, F] = mesh_filter_vertices(dualquat.vertices, dualquat.faces, filter);
% clipping_plane = V(:, 3) > 0;
% V(clipping_plane, :) = NaN(nnz(clipping_plane), 3);
trimesh(F, V(:, 1), V(:, 2), V(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 1.0, 'SpecularStrength', 0.4);
view(l1); camlight;
view(l2); camlight;
lighting gouraud;
view(vp);
zoom(1.2);
axis equal;
grid off;
axis off;
pause;
print('figs/index_dualquat.png', '-dpng', '-r300');

% Use implicit skinning
implicit = skin_implicit(mesh, transforms, 1, 'hrbf');
[V, F] = mesh_filter_vertices(implicit.vertices, implicit.faces, filter);
% clipping_plane = V(:, 3) > 2;
% V(clipping_plane, :) = NaN(nnz(clipping_plane), 3);
trimesh(F, V(:, 1), V(:, 2), V(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 1.0, 'SpecularStrength', 0.4);
view(l1); camlight;
view(l2); camlight;
lighting gouraud;
view(vp);
zoom(1.2);
axis equal;
grid off;
axis off;
pause;
print('figs/index_implicit.png', '-dpng', '-r300');

implicit = skin_implicit(mesh, transforms, 3, 'hrbf');
[V, F] = mesh_filter_vertices(implicit.vertices, implicit.faces, filter);
% clipping_plane = V(:, 3) > 2;
% V(clipping_plane, :) = NaN(nnz(clipping_plane), 3);
trimesh(F, V(:, 1), V(:, 2), V(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 1.0, 'SpecularStrength', 0.4);
view(l1); camlight;
view(l2); camlight;
lighting gouraud;
view(vp);
zoom(1.2);
axis equal;
grid off;
axis off;
pause;
print('figs/index_implicit_thresh3.png', '-dpng', '-r300');

reference = {};
[reference.vertices, reference.faces] = load_mesh('mesh/3.obj');
load('mesh/centers3.mat');
centers = zeros(length(centers3), 3);
for i = 1 : length(centers3)
    centers(i, :) = [centers3{i}(2), centers3{i}(3), centers3{i}(1)];
end
transforms = bone_transforms(mesh, centers);

implicit = skin_implicit(mesh, transforms, 1, 'spheres');
trimesh(implicit.faces, implicit.vertices(:, 1), implicit.vertices(:, 2), implicit.vertices(:, 3), 'EdgeColor', [0.4, 0.25, 0.0], 'EdgeAlpha', 0.0, 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 1.0, 'SpecularStrength', 0.4);
view([-25, 25]);
camlight;
view([0, -65]);
camlight;
axis equal;
grid off;
lighting gouraud;
axis off;
zoom(2);
pause;
print('figs/palm_implicit.png', '-dpng', '-r300');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;

% Load mesh
load('mesh/neutral.mat');
l1 = [180, 40];
l2 = [180, -30];
vp = [180, 10];
z = 1.4;
colors = hsv(18);

% Show spheres
clf;
trimesh(mesh.faces, mesh.vertices(:, 1), mesh.vertices(:, 2), mesh.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.5);
hold on;
for i = 1 : 27
    [X, Y, Z] = sphere(100);
    X = X * mesh.spheres{i}.radius + mesh.spheres{i}.center(1);
    Y = Y * mesh.spheres{i}.radius + mesh.spheres{i}.center(2);
    Z = Z * mesh.spheres{i}.radius + mesh.spheres{i}.center(3);
    surf(X, Y, Z, 'FaceColor', [0.2, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    text(mesh.spheres{i}.center(1), mesh.spheres{i}.center(2), mesh.spheres{i}.center(3), num2str(i), 'Color', [0, 0, 0], 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end
for i = 1 : length(mesh.blocks)
    for j = 1 : length(mesh.blocks{i})
        a = mesh.spheres{mesh.blocks{i}(j)}.center;
        b = mesh.spheres{mesh.blocks{i}(mod(j, length(mesh.blocks{i})) + 1)}.center;
        plot3([a(1), b(1)], [a(2), b(2)], [a(3), b(3)], 'Color', [0.5, 0.5, 0.5], 'LineWidth', 3);
    end
end
hold off;
view(l1); camlight;
view(l2); camlight;
view(vp);
axis equal;
grid off;
lighting gouraud;
axis off;
zoom(z);
pause;
print('figs/hand_spheres.png', '-dpng', '-r300');

% Show bones
clf;
axes = bone_axes(mesh.spheres);
trimesh(mesh.faces, mesh.vertices(:, 1), mesh.vertices(:, 2), mesh.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.5);
hold on;
for i = 1 : length(mesh.bones)
    [X, Y, Z] = sphere(100);
    X = X * 0.75 + mesh.bones{i}.center(1);
    Y = Y * 0.75 + mesh.bones{i}.center(2);
    Z = Z * 0.75 + mesh.bones{i}.center(3);
    surf(X, Y, Z, 'EdgeColor', 'none', 'FaceColor', colors(i, :), 'FaceAlpha', 0.5);
    if mesh.bones{i}.parent ~= 0
        a = mesh.bones{i}.center;
        b = mesh.bones{mesh.bones{i}.parent}.center;
        plot3([a(1), b(1)], [a(2), b(2)], [a(3), b(3)], 'Color', [0.5, 0.5, 0.5], 'LineWidth', 3);
    end
    quiver3(repmat(axes{i}(1, 4), 1, 3), repmat(axes{i}(2, 4), 1, 3), repmat(axes{i}(3, 4), 1, 3), axes{i}(1, 1 : 3), axes{i}(2, 1 : 3), axes{i}(3, 1 : 3), 'Color', colors(i, :), 'LineWidth', 2);
    text(mesh.bones{i}.center(1), mesh.bones{i}.center(2), mesh.bones{i}.center(3), num2str(i), 'Color', [0, 0, 0], 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end
hold off;
view(l1); camlight;
view(l2); camlight;
view(vp);
axis equal;
grid off;
lighting gouraud;
axis off;
zoom(z);
pause;
print('figs/hand_bones.png', '-dpng', '-r300');

% Show bone weights
clf;
patch('Faces', mesh.faces, 'Vertices', mesh.vertices, 'FaceVertexCData', mesh.weights * colors, ...
    'EdgeColor', 'none', 'FaceColor', 'interp', 'CDataMapping', 'scaled')
view(l1); camlight;
view(l2); camlight;
view(vp);
axis equal;
grid off;
lighting gouraud;
axis off;
zoom(z);
pause;
print('figs/hand_bone_weights.png', '-dpng', '-r300');

% Add spacing between bones
offset = 3;
delta = zeros(18, 3);
for i = 1 : 18
    delta(i, :) = offset * axes{i}(1 : 3, 1)';
    if mesh.bones{i}.parent > 0
        delta(i, :) = delta(i, :) + delta(mesh.bones{i}.parent, :);
    end
end

% Show bone partitioning
clf;
for i = 1 : 18
    if i == 2
        hold on;
    end
    patch('Faces', mesh.bones{i}.faces, 'Vertices', mesh.bones{i}.vertices + repmat(delta(i, :), size(mesh.bones{i}.vertices, 1), 1), 'FaceVertexCData', repmat(colors(i, :), size(mesh.bones{i}.vertices, 1), 1), ...
        'EdgeColor', 'none', 'FaceColor', 'interp', 'CDataMapping', 'scaled')
end
hold off;
view(l1); camlight;
view(l2); camlight;
view(vp);
axis equal;
grid off;
lighting gouraud;
axis off;
zoom(z);
pause;
print('figs/hand_bone_partition.png', '-dpng', '-r300');

% Show bone HRBF reconstruction
clf;
for i = 1 : 18
    if i == 2
        hold on;
    end
    [X, Y, Z] = meshgrid( ...
        min(mesh.bones{i}.vertices(:, 1)) - 2 : 0.1 : max(mesh.bones{i}.vertices(:, 1)) + 3, ...
        min(mesh.bones{i}.vertices(:, 2)) - 2 : 0.1 : max(mesh.bones{i}.vertices(:, 2)) + 2, ...
        min(mesh.bones{i}.vertices(:, 3)) - 2 : 0.1 : max(mesh.bones{i}.vertices(:, 3)) + 2 ...
    );
    V = hrbf_apply(mesh.bones{i}.hrbf.centers, mesh.bones{i}.hrbf.coefficients, [X(:), Y(:), Z(:)]);
    V = reshape(V, size(X));
    X = X + delta(i, 1);
    Y = Y + delta(i, 2);
    Z = Z + delta(i, 3);
    p = patch(isosurface(X, Y, Z, V, 0));
    isonormals(X, Y, Z, V, p);
    p.FaceColor = colors(i, :);
    p.EdgeColor = 'none';
end
hold off;
view(l1); camlight;
view(l2); camlight;
view(vp);
axis equal;
grid off;
lighting gouraud;
axis off;
zoom(z);
pause;
print('figs/hand_bone_hrbf.png', '-dpng', '-r300');

% Show bone sphere reconstruction
clf;
for i = 1 : 18
    if i == 2
        hold on;
    end
    [X, Y, Z] = meshgrid( ...
        min(mesh.bones{i}.vertices(:, 1)) - 2 : 0.1 : max(mesh.bones{i}.vertices(:, 1)) + 2, ...
        min(mesh.bones{i}.vertices(:, 2)) - 2 : 0.1 : max(mesh.bones{i}.vertices(:, 2)) + 2, ...
        min(mesh.bones{i}.vertices(:, 3)) - 2 : 0.1 : max(mesh.bones{i}.vertices(:, 3)) + 2 ...
    );
    V = sphere_apply(mesh.spheres, mesh.bones{i}.blocks, [X(:), Y(:), Z(:)]);
    V = reshape(V, size(X));
    X = X + delta(i, 1);
    Y = Y + delta(i, 2);
    Z = Z + delta(i, 3);
    p = patch(isosurface(X, Y, Z, V, 0));
    isonormals(X, Y, Z, V, p);
    p.FaceColor = colors(i, :);
    p.EdgeColor = 'none';
end
hold off;
view(l1); camlight;
view(l2); camlight;
view(vp);
axis equal;
grid off;
lighting gouraud;
axis off;
zoom(z);
pause;
print('figs/hand_bone_sphere.png', '-dpng', '-r300');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;

% Load mesh
load('mesh/neutral.mat');
l1 = [0, 60];
l2 = [0, -60];
vp = [0, 10];
z = 1.4;
colors = hsv(18);

% Show bone weights
clf;
patch('Faces', mesh.faces, 'Vertices', mesh.vertices, 'FaceVertexCData', ...
    mesh.implicit.blend_weight * [1.0, 0.7, 0.5] + (1 - mesh.implicit.blend_weight) * [0.8, 0.8, 0.8], ...
    'EdgeColor', 'none', 'FaceColor', 'interp', 'CDataMapping', 'scaled')
view(l1); camlight;
view(l2); camlight;
view(vp);
axis equal;
grid off;
lighting gouraud;
axis off;
zoom(z);
pause;
print('figs/hand_bone_blend_weights.png', '-dpng', '-r300');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;

% Load mesh
load('mesh/neutral.mat');
l1 = [0, 30];
l2 = [180, 30];
vp = [90, 90];
z = 1.1;
threshold = 4;

% Keep only the index
bones = mesh.bones(6 : 8);
bones{1}.parent = 0;
bones{2}.parent = 1;
bones{3}.parent = 2;
transforms = {eye(4), [], []};
transforms{2} = matrix_rotation(pi / 3, [0, 0, 1], bones{2}.center);
transforms{3} = transforms{2} * matrix_rotation(pi / 2, [0, 0, 1], bones{3}.center);
[X, Y, Z] = meshgrid( ...
    -9 : 0.1 : 0, ...
    -3 : 0.1 : 5, ...
    -1 : 0.1 : 5 ...
);

% Show HRBF reconstruction with union
clf;
V1 = hrbf_apply(bones{1}.hrbf.centers, bones{1}.hrbf.coefficients, matrix_apply(inv(transforms{1}), [X(:), Y(:), Z(:)]));
V2 = hrbf_apply(bones{2}.hrbf.centers, bones{2}.hrbf.coefficients, matrix_apply(inv(transforms{2}), [X(:), Y(:), Z(:)]));
V3 = hrbf_apply(bones{3}.hrbf.centers, bones{3}.hrbf.coefficients, matrix_apply(inv(transforms{3}), [X(:), Y(:), Z(:)]));
V1 = sdf_reparam(V1, threshold);
V2 = sdf_reparam(V2, threshold);
V3 = sdf_reparam(V3, threshold);
V12 = sdf_union(V1, V2);
V23 = sdf_union(V2, V3);
V = sdf_union(V12, V23);
V = reshape(V, size(X));
p = patch(isosurface(X, Y, Z, V, 0.5));
isonormals(X, Y, Z, V, p);
p.FaceColor = [0.5, 0.5, 0.8];
p.FaceAlpha = 0.5;
p.EdgeColor = 'none';
for s = [0.4, 0.3, 0.2, 0.1]
    p = patch(isosurface(X, Y, Z, V, s));
    isonormals(X, Y, Z, V, p);
    p.FaceColor = [0.5, 0.5, 0.6];
    p.FaceAlpha = 0.1;
    p.EdgeColor = 'none';
end
view(l1); camlight;
view(l2); camlight;
lighting gouraud;
view(vp);
zoom(z);
axis equal;
grid off;
axis off;
pause;
print('figs/index_sdf_hrbf_union.png', '-dpng', '-r300');

% Show sphere reconstruction with union
clf;
V1 = sphere_apply(mesh.spheres, bones{1}.blocks, matrix_apply(inv(transforms{1}), [X(:), Y(:), Z(:)]));
V2 = sphere_apply(mesh.spheres, bones{2}.blocks, matrix_apply(inv(transforms{2}), [X(:), Y(:), Z(:)]));
V3 = sphere_apply(mesh.spheres, bones{3}.blocks, matrix_apply(inv(transforms{3}), [X(:), Y(:), Z(:)]));
V1 = sdf_reparam(V1, threshold);
V2 = sdf_reparam(V2, threshold);
V3 = sdf_reparam(V3, threshold);
V12 = sdf_union(V1, V2);
V23 = sdf_union(V2, V3);
V = sdf_union(V12, V23);
V = reshape(V, size(X));
p = patch(isosurface(X, Y, Z, V, 0.5));
isonormals(X, Y, Z, V, p);
p.FaceColor = [0.5, 0.5, 0.8];
p.FaceAlpha = 0.5;
p.EdgeColor = 'none';
for s = [0.4, 0.3, 0.2, 0.1]
    p = patch(isosurface(X, Y, Z, V, s));
    isonormals(X, Y, Z, V, p);
    p.FaceColor = [0.5, 0.5, 0.6];
    p.FaceAlpha = 0.1;
    p.EdgeColor = 'none';
end
view(l1); camlight;
view(l2); camlight;
lighting gouraud;
view(vp);
zoom(z);
axis equal;
grid off;
axis off;
pause;
print('figs/index_sdf_sphere_union.png', '-dpng', '-r300');

% Show HRBF reconstruction with blending
clf;
alpha_to_theta = {};
alpha_to_theta.alpha0 = pi / 8;
alpha_to_theta.alpha1 = pi / 2;
alpha_to_theta.alpha2 = pi / 2;
alpha_to_theta.theta0 = pi / 4;
alpha_to_theta.theta1 = pi / 10;
alpha_to_theta.theta2 = pi / 10;
alpha_to_theta.omega0 = 0.5;
alpha_to_theta.omega1 = 1;
[V1, G1] = hrbf_apply(bones{1}.hrbf.centers, bones{1}.hrbf.coefficients, matrix_apply(inv(transforms{1}), [X(:), Y(:), Z(:)]));
[V2, G2] = hrbf_apply(bones{2}.hrbf.centers, bones{2}.hrbf.coefficients, matrix_apply(inv(transforms{2}), [X(:), Y(:), Z(:)]));
[V3, G3] = hrbf_apply(bones{3}.hrbf.centers, bones{3}.hrbf.coefficients, matrix_apply(inv(transforms{3}), [X(:), Y(:), Z(:)]));
V1 = sdf_reparam(V1, threshold);
V2 = sdf_reparam(V2, threshold);
V3 = sdf_reparam(V3, threshold);
G1 = -normalizerow(G1);
G2 = -normalizerow(G2);
G3 = -normalizerow(G3);
% alpha = acos(dot(G1, G2, 2));
% theta = sdf_blend_theta_curve(alpha, alpha_to_theta);
% V12 = sdf_blend(V1, V2, theta);
% alpha = acos(dot(G2, G3, 2));
% theta = sdf_blend_theta_curve(alpha, alpha_to_theta);
% V23 = sdf_blend(V2, V3, theta);
V12b = sdf_blend(V1, V2, sdf_blend_theta_curve(pi / 3, alpha_to_theta));
V23b = sdf_blend(V2, V3, sdf_blend_theta_curve(pi / 2, alpha_to_theta));
V12u = sdf_union(V1, V2);
V23u = sdf_union(V2, V3);
% TODO blending weights
V12 = V12b;
V23 = V23b;
V = sdf_union(V12, V23);
V = reshape(V, size(X));
p = patch(isosurface(X, Y, Z, V, 0.5));
isonormals(X, Y, Z, V, p);
p.FaceColor = [0.5, 0.5, 0.8];
p.FaceAlpha = 0.5;
p.EdgeColor = 'none';
for s = [0.4, 0.3, 0.2, 0.1]
    p = patch(isosurface(X, Y, Z, V, s));
    isonormals(X, Y, Z, V, p);
    p.FaceColor = [0.5, 0.5, 0.6];
    p.FaceAlpha = 0.1;
    p.EdgeColor = 'none';
end
view(l1); camlight;
view(l2); camlight;
lighting gouraud;
view(vp);
zoom(z);
axis equal;
grid off;
axis off;
pause;
print('figs/index_sdf_hrbf_blend.png', '-dpng', '-r300');

% Show sphere reconstruction with blending
clf;
[V1, G1] = sphere_apply(mesh.spheres, bones{1}.blocks, matrix_apply(inv(transforms{1}), [X(:), Y(:), Z(:)]));
[V2, G2] = sphere_apply(mesh.spheres, bones{2}.blocks, matrix_apply(inv(transforms{2}), [X(:), Y(:), Z(:)]));
[V3, G3] = sphere_apply(mesh.spheres, bones{3}.blocks, matrix_apply(inv(transforms{3}), [X(:), Y(:), Z(:)]));
V1 = sdf_reparam(V1, threshold);
V2 = sdf_reparam(V2, threshold);
V3 = sdf_reparam(V3, threshold);
G1 = -normalizerow(G1);
G2 = -normalizerow(G2);
G3 = -normalizerow(G3);
% alpha = acos(dot(G1, G2, 2));
% theta = sdf_blend_theta_curve(alpha, alpha_to_theta);
% V12 = sdf_blend(V1, V2, theta);
% alpha = acos(dot(G2, G3, 2));
% theta = sdf_blend_theta_curve(alpha, alpha_to_theta);
% V23 = sdf_blend(V2, V3, theta);
V12b = sdf_blend(V1, V2, sdf_blend_theta_curve(pi / 3, alpha_to_theta));
V23b = sdf_blend(V2, V3, sdf_blend_theta_curve(pi / 2, alpha_to_theta));
V12u = sdf_union(V1, V2);
V23u = sdf_union(V2, V3);
% TODO blending weights
V12 = V12b;
V23 = V23b;
V = sdf_union(V12, V23);
V = real(reshape(V, size(X)));
p = patch(isosurface(X, Y, Z, V, 0.5));
isonormals(X, Y, Z, V, p);
p.FaceColor = [0.5, 0.5, 0.8];
p.FaceAlpha = 0.5;
p.EdgeColor = 'none';
for s = [0.4, 0.3, 0.2, 0.1]
    p = patch(isosurface(X, Y, Z, V, s));
    isonormals(X, Y, Z, V, p);
    p.FaceColor = [0.5, 0.5, 0.6];
    p.FaceAlpha = 0.1;
    p.EdgeColor = 'none';
end
view(l1); camlight;
view(l2); camlight;
lighting gouraud;
view(vp);
zoom(z);
axis equal;
grid off;
axis off;
pause;
print('figs/index_sdf_sphere_blend.png', '-dpng', '-r300');
