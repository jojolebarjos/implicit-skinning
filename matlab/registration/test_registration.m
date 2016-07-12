
% Register gptoolbox
addpath(genpath('../external'));

% Load neutral mesh and point cloud
load('../mesh/neutral.mat');
points = {};
[points.vertices, points.faces] = load_mesh('anastasia.obj');
points.normals = per_vertex_normals(points.vertices, points.faces);
transforms = cell(1, 18);
for i = 1 : 18
    transforms{i} = eye(4);
end
axes = compute_bone_axes(mesh.spheres);

% Compute initial guess using while mesh and cloud
transform = inv(initial_guess(points.vertices)) * initial_guess(mesh.vertices);
faces = mesh.faces;
vertices = apply_matrix(transform, mesh.vertices);
normals = -per_vertex_normals(vertices, faces);
trimesh(faces, vertices(:, 1), vertices(:, 2), vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.9, 0.4], 'FaceAlpha', 0.5);
hold on;
trimesh(points.faces, points.vertices(:, 1), points.vertices(:, 2), points.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.5);
quiver3(vertices(:, 1), vertices(:, 2), vertices(:, 3), normals(:, 1), normals(:, 2), normals(:, 3), 'Color', [0.4, 0.9, 0.4]);
quiver3(points.vertices(:, 1), points.vertices(:, 2), points.vertices(:, 3), points.normals(:, 1), points.normals(:, 2), points.normals(:, 3), 'Color', [0.8, 0.8, 0.8]);
hold off;
view([-90, 0]);
camlight;
view([90, 0]);
camlight;
axis equal;
grid off;
lighting gouraud;
axis off;
title('Initial guess');
pause;

% Keep only the palm and find rigid transform iteratively
keep = ismember(mesh.assignments, [2, 3, 6, 9, 12, 15, 16]);
[vertices, faces] = filter_vertices(vertices, faces, keep);
normals = normals(keep, :);
pairs = compute_correspondences(vertices, normals, points.vertices, points.normals);
for i = 1 : 20
    delta = compute_transformation(vertices, points.vertices, points.normals, pairs);
    transform = delta * transform;
    vertices = apply_matrix(delta, vertices);
    normals = apply_matrix(delta, normals, 0);
    pairs = compute_correspondences(vertices, normals, points.vertices, points.normals);
    v = get(gca, 'view');
    trimesh(faces, vertices(:, 1), vertices(:, 2), vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.9, 0.4], 'FaceAlpha', 0.1);
    hold on;
    trimesh(points.faces, points.vertices(:, 1), points.vertices(:, 2), points.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.1);
    plot3( ...
        [vertices(pairs(:, 1), 1), points.vertices(pairs(:, 2), 1)]', ...
        [vertices(pairs(:, 1), 2), points.vertices(pairs(:, 2), 2)]', ...
        [vertices(pairs(:, 1), 3), points.vertices(pairs(:, 2), 3)]', ...
    'Color', 'red');
    hold off;
    view([-90, 0]);
    camlight;
    view([90, 0]);
    camlight;
    axis equal;
    grid off;
    lighting gouraud;
    axis off;
    title(['After ', num2str(i), ' rigid transformation']);
    set(gca, 'view', v);
    pause(0.01);
end
transforms{2} = transform;
% TODO improve proto-phalange?
transforms{15} = transform;
%pause;

% Fit finger roots
bones = [6, 9, 12, 16];
for f = 1 : 4
    transform = transforms{2};
    keep = mesh.assignments == bones(f);
    [vertices, faces] = filter_vertices(mesh.vertices, mesh.faces, keep);
    normals = mesh.normals(keep, :);
    vertices = apply_matrix(transform, vertices);
    normals = apply_matrix(transform, normals, 0);
    pairs = compute_correspondences(vertices, normals, points.vertices, points.normals, 0.1, 0.2);
    for i = 1 : 10
        delta = compute_transformation(vertices, points.vertices, points.normals, pairs);
        delta = constraint_transformation(delta, apply_matrix(transforms{2}, mesh.bones{bones(f)}.center));
        transform = delta * transform;
        vertices = apply_matrix(delta, vertices);
        normals = apply_matrix(delta, normals, 0);
        pairs = compute_correspondences(vertices, normals, points.vertices, points.normals, 0.3, 0.3);
        v = get(gca, 'view');
        trimesh(faces, vertices(:, 1), vertices(:, 2), vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.9, 0.4], 'FaceAlpha', 0.1);
        hold on;
        trimesh(points.faces, points.vertices(:, 1), points.vertices(:, 2), points.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.1);
        plot3( ...
            [vertices(pairs(:, 1), 1), points.vertices(pairs(:, 2), 1)]', ...
            [vertices(pairs(:, 1), 2), points.vertices(pairs(:, 2), 2)]', ...
            [vertices(pairs(:, 1), 3), points.vertices(pairs(:, 2), 3)]', ...
        'Color', 'red');
        hold off;
        view([-90, 0]);
        camlight;
        view([90, 0]);
        camlight;
        axis equal;
        grid off;
        lighting gouraud;
        axis off;
        title(['After ', num2str(i), ' rigid transformation']);
        set(gca, 'view', v);
        pause(0.01);
    end
    transforms{bones(f)} = transform;
    %pause;
end

% Improve palm origin
% transform = transforms{2};
% keep = ismember(mesh.assignments, [2, 3, 15]);
% [vertices, faces] = filter_vertices(mesh.vertices, mesh.faces, keep);
% normals = mesh.normals(keep, :);
% vertices = apply_matrix(transform, vertices);
% normals = apply_matrix(transform, normals, 0);
% pairs = compute_correspondences(vertices, normals, points.vertices, points.normals);
% for i = 1 : 10
%     delta = compute_transformation(vertices, points.vertices, points.normals, pairs);
%     %delta = constraint_transformation(delta);
%     transform = delta * transform;
%     vertices = apply_matrix(delta, vertices);
%     normals = apply_matrix(delta, normals, 0);
%     pairs = compute_correspondences(vertices, normals, points.vertices, points.normals);
%     v = get(gca, 'view');
%     trimesh(faces, vertices(:, 1), vertices(:, 2), vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.9, 0.4], 'FaceAlpha', 0.1);
%     hold on;
%     trimesh(points.faces, points.vertices(:, 1), points.vertices(:, 2), points.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.1);
%     plot3( ...
%         [vertices(pairs(:, 1), 1), points.vertices(pairs(:, 2), 1)]', ...
%         [vertices(pairs(:, 1), 2), points.vertices(pairs(:, 2), 2)]', ...
%         [vertices(pairs(:, 1), 3), points.vertices(pairs(:, 2), 3)]', ...
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
%     pause(0.01);
% end
% transforms{2} = transform;
% % TODO improve proto-phalange?
% transforms{15} = transform;
% pause;

% Estimate wrist rotation
transform = transforms{2};
keep = mesh.assignments == 1;
[vertices, faces] = filter_vertices(mesh.vertices, mesh.faces, keep);
normals = mesh.normals(keep, :);
vertices = apply_matrix(transform, vertices);
normals = apply_matrix(transform, normals, 0);
pairs = compute_correspondences(vertices, normals, points.vertices, points.normals, 0.1, 0.2);
for i = 1 : 10
    delta = compute_transformation(vertices, points.vertices, points.normals, pairs);
    delta = constraint_transformation(delta, apply_matrix(transforms{2}, mesh.bones{2}.center));
    transform = delta * transform;
    vertices = apply_matrix(delta, vertices);
    normals = apply_matrix(delta, normals, 0);
    pairs = compute_correspondences(vertices, normals, points.vertices, points.normals, 0.3, 0.3);
    v = get(gca, 'view');
    trimesh(faces, vertices(:, 1), vertices(:, 2), vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.9, 0.4], 'FaceAlpha', 0.1);
    hold on;
    trimesh(points.faces, points.vertices(:, 1), points.vertices(:, 2), points.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.1);
    plot3( ...
        [vertices(pairs(:, 1), 1), points.vertices(pairs(:, 2), 1)]', ...
        [vertices(pairs(:, 1), 2), points.vertices(pairs(:, 2), 2)]', ...
        [vertices(pairs(:, 1), 3), points.vertices(pairs(:, 2), 3)]', ...
    'Color', 'red');
    hold off;
    view([-90, 0]);
    camlight;
    view([90, 0]);
    camlight;
    axis equal;
    grid off;
    lighting gouraud;
    axis off;
    title(['After ', num2str(i), ' rigid transformation']);
    set(gca, 'view', v);
    pause(0.01);
end
transforms{1} = transform;
%pause;

% Fit thumb root
transform = transforms{2};
keep = ismember(mesh.assignments, [3, 4, 5]);
[vertices, faces] = filter_vertices(mesh.vertices, mesh.faces, keep);
normals = mesh.normals(keep, :);
vertices = apply_matrix(transform, vertices);
normals = apply_matrix(transform, normals, 0);
pairs = compute_correspondences(vertices, normals, points.vertices, points.normals, 0.1, 0.2);
for i = 1 : 10
    delta = compute_transformation(vertices, points.vertices, points.normals, pairs);
    delta = constraint_transformation(delta, apply_matrix(transforms{2}, mesh.bones{3}.center), apply_matrix(transforms{2}, axes{2}(1 : 3, 2)', 0));
    transform = delta * transform;
    vertices = apply_matrix(delta, vertices);
    normals = apply_matrix(delta, normals, 0);
    pairs = compute_correspondences(vertices, normals, points.vertices, points.normals, 0.3, 0.3);
    v = get(gca, 'view');
    trimesh(faces, vertices(:, 1), vertices(:, 2), vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.9, 0.4], 'FaceAlpha', 0.1);
    hold on;
    trimesh(points.faces, points.vertices(:, 1), points.vertices(:, 2), points.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.1);
    plot3( ...
        [vertices(pairs(:, 1), 1), points.vertices(pairs(:, 2), 1)]', ...
        [vertices(pairs(:, 1), 2), points.vertices(pairs(:, 2), 2)]', ...
        [vertices(pairs(:, 1), 3), points.vertices(pairs(:, 2), 3)]', ...
    'Color', 'red');
    hold off;
    view([-90, 0]);
    camlight;
    view([90, 0]);
    camlight;
    axis equal;
    grid off;
    lighting gouraud;
    axis off;
    title(['After ', num2str(i), ' rigid transformation']);
    set(gca, 'view', v);
    pause(0.01);
end
transforms{3} = transform;
transforms{4} = transform;
%pause;

% Use delineation to estime finger lengths and pose
bones = [6, 9, 12, 16];
for f = bones
    
end

% TODO recompute skeleton

% Show registered bones
colors = hsv(18);
trimesh(points.faces, points.vertices(:, 1), points.vertices(:, 2), points.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.1);
hold on;
for i = 1 : 18
    if transforms{i}(1, 1) == 1
        continue;
    end
    keep = mesh.assignments == i;
    [vertices, faces] = filter_vertices(mesh.vertices, mesh.faces, keep);
    normals = mesh.normals(keep, :);
    vertices = apply_matrix(transforms{i}, vertices);
    normals = apply_matrix(transforms{i}, normals, 0);
    trimesh(faces, vertices(:, 1), vertices(:, 2), vertices(:, 3), 'EdgeColor', colors(i, :), 'FaceColor', 'none');
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
title(['After ', num2str(i), ' rigid transformation']);
set(gca, 'view', v);
