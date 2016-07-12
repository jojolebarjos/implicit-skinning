
% Register gptoolbox
addpath(genpath('../external'));

% Load cleaned mesh and point cloud
cleaned = {};
[cleaned.vertices, cleaned.faces] = load_mesh('anastasia_clean.obj');
points = {};
[points.vertices, points.faces] = load_mesh('anastasia_cloud.obj');
points.normals = per_vertex_normals(points.vertices, points.faces);

% Estimate pose
% transform = inv(initial_guess(points.vertices)) * [-0.5,0,0,-1;0,0.5,0,0;0,0,0.5,0;0,0,0,1] * initial_guess(cleaned.vertices);
% cleaned.vertices = apply_matrix(transform, cleaned.vertices);
% cleaned.normals = per_vertex_normals(cleaned.vertices, cleaned.faces);
% trimesh(cleaned.faces, cleaned.vertices(:, 1), cleaned.vertices(:, 2), cleaned.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.9, 0.4], 'FaceAlpha', 0.5);
% hold on;
% trimesh(points.faces, points.vertices(:, 1), points.vertices(:, 2), points.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.5);
% quiver3(cleaned.vertices(:, 1), cleaned.vertices(:, 2), cleaned.vertices(:, 3), cleaned.normals(:, 1), cleaned.normals(:, 2), cleaned.normals(:, 3), 'Color', [0.4, 0.9, 0.4]);
% quiver3(points.vertices(:, 1), points.vertices(:, 2), points.vertices(:, 3), points.normals(:, 1), points.normals(:, 2), points.normals(:, 3), 'Color', [0.8, 0.8, 0.8]);
% hold off;
% view([-90, 0]);
% camlight;
% view([90, 0]);
% camlight;
% axis equal;
% grid off;
% lighting gouraud;
% axis off;
% title('Initial guess');
% pause;

% Find rigid transform iteratively
% pairs = compute_correspondences(cleaned.vertices, cleaned.normals, points.vertices, points.normals);
% for i = 1 : 10
%     delta = compute_transformation(cleaned.vertices, points.vertices, points.normals, pairs);
%     transform = delta * transform;
%     cleaned.vertices = apply_matrix(delta, cleaned.vertices);
%     cleaned.normals = apply_matrix(delta, cleaned.normals, 0);
%     pairs = compute_correspondences(cleaned.vertices, cleaned.normals, points.vertices, points.normals);
%     v = get(gca, 'view');
%     trimesh(cleaned.faces, cleaned.vertices(:, 1), cleaned.vertices(:, 2), cleaned.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.9, 0.4], 'FaceAlpha', 0.1);
%     hold on;
%     trimesh(points.faces, points.vertices(:, 1), points.vertices(:, 2), points.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.1);
%     plot3( ...
%         [cleaned.vertices(pairs(:, 1), 1), points.vertices(pairs(:, 2), 1)]', ...
%         [cleaned.vertices(pairs(:, 1), 2), points.vertices(pairs(:, 2), 2)]', ...
%         [cleaned.vertices(pairs(:, 1), 3), points.vertices(pairs(:, 2), 3)]', ...
%     'Color', 'red');
%     hold off;
%     view([-90, 0]);
%     camlight;
%     view([90, 0]);
%     camlight;
%     axis equal;
%     grid off;
%     lighting gouraud;
%     axis off;
%     title(['After ', num2str(i), ' rigid transformation']);
%     set(gca, 'view', v);
%     pause;
% end
transform = [ ...
    0.1607   -0.2572    0.0193   -0.8139; ...
    0.2539    0.1538   -0.0650   -0.4930; ...
    0.0453    0.0505    0.2962    7.2152; ...
         0         0         0    1.0000];

% Load and transform sphere locations
load('anastasia_clean.mat');
% Note: spheres do NOT have the same topology as before...
% 1 to 20 do not change
% 21 -> 25
% 22 -> 26
% 23 -> 28
% 24 -> 30
% 25 -> 29
% 26 -> 31
% 27 -> 32
map = [1 : 20, 25, 26, 28, 30, 29, 31, 32];
centers = centers(map, :);
radii = radii(map);
centers = apply_matrix(transform, centers);
radii = radii * norm(transform(1 : 3, 1));

% Flip everything, since neutral mesh is right-handed
axe = princomp(points.vertices);
normal = axe(2, :);
reflect = eye(3) - 2 * (normal' * normal);
points.vertices = apply_matrix(reflect, points.vertices);
points.normals = apply_matrix(reflect, points.normals, 0);
centers = apply_matrix(reflect, centers);

% Show spheres
v = get(gca, 'view');
trimesh(points.faces, points.vertices(:, 1), points.vertices(:, 2), points.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.5);
hold on;
for i = 1 : length(radii)
    [X, Y, Z] = sphere;
    X = X * radii(i) + centers(i, 1);
    Y = Y * radii(i) + centers(i, 2);
    Z = Z * radii(i) + centers(i, 3);
    surf(X, Y, Z, 'FaceColor', [0.2, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    text(centers(i, 1), centers(i, 2), centers(i, 3), num2str(i), 'Color', [0, 0, 0], 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end
hold off;
view([-90, 0]);
camlight;
view([90, 0]);
camlight;
axis equal;
grid off;
lighting gouraud;
axis off;
set(gca, 'view', v);
pause;

% Show bones
v = get(gca, 'view');
colors = hsv(18);
axes = compute_bone_axes(centers);
trimesh(points.faces, points.vertices(:, 1), points.vertices(:, 2), points.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.5);
hold on;
for i = 1 : 18
    [X, Y, Z] = sphere;
    X = X * 0.05 + axes{i}(1, 4);
    Y = Y * 0.05 + axes{i}(2, 4);
    Z = Z * 0.05 + axes{i}(3, 4);
    surf(X, Y, Z, 'EdgeColor', 'none', 'FaceColor', colors(i, :), 'FaceAlpha', 0.5);
    quiver3(repmat(axes{i}(1, 4), 1, 3), repmat(axes{i}(2, 4), 1, 3), repmat(axes{i}(3, 4), 1, 3), 0.25 * axes{i}(1, 1 : 3), 0.25 * axes{i}(2, 1 : 3), 0.25 * axes{i}(3, 1 : 3), 'Color', colors(i, :), 'LineWidth', 2);
    text(axes{i}(1, 4), axes{i}(2, 4), axes{i}(3, 4), num2str(i), 'Color', [0, 0, 0], 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end
hold off;
view([-90, 0]);
camlight;
view([90, 0]);
camlight;
axis equal;
grid off;
lighting gouraud;
axis off;
set(gca, 'view', v);
pause;

% Use linear skinning to approximate point cloud
load('../final/neutral.mat');
transforms = compute_transfer_transforms(mesh, centers, radii);
mesh = skin_linear(mesh, transforms);
%mesh.vertices = apply_matrix(transforms{3}, mesh.vertices);
trimesh(mesh.faces, mesh.vertices(:, 1), mesh.vertices(:, 2), mesh.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.9, 0.4], 'FaceAlpha', 0.5);
hold on;
trimesh(points.faces, points.vertices(:, 1), points.vertices(:, 2), points.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.5);
hold off;
view([-90, 0]);
camlight;
view([90, 0]);
camlight;
axis equal;
grid off;
lighting gouraud;
axis off;
