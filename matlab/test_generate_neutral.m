
% Register gptoolbox
addpath(genpath('external'));

% Load mesh
mesh = {};
[mesh.vertices, mesh.faces, ~, ~, N, NF] = readOBJ('mesh/5.obj');
M = zeros(size(mesh.vertices, 1), 1);
for f = 1 : size(mesh.faces, 1)
    M(mesh.faces(f, 1)) = NF(f, 1);
    M(mesh.faces(f, 2)) = NF(f, 2);
    M(mesh.faces(f, 3)) = NF(f, 3);
end
mesh.normals = N(M, :);

% Load and rescale additional poses
%load('mesh/centers1.mat');
load('mesh/centers2.mat');
load('mesh/centers3.mat');
%load('mesh/centers_spock.mat');
load('mesh/centers_ok.mat');
load('mesh/centers_scissors.mat');
spheres = {centers2, centers3, centers_ok, centers_scissors};
poses = cell(size(spheres));
for i = 1 : length(spheres)
    centers = zeros(30, 3);
    for j = 1 : 30
        centers(j, 1) = spheres{i}{j}(2);
        centers(j, 2) = spheres{i}{j}(3);
        centers(j, 3) = spheres{i}{j}(1);
    end
    % TODO rescale poses
    poses{i} = centers;
end
save('mesh/poses.mat', 'poses');
disp('poses.mat saved');

% Show mesh
trimesh(mesh.faces, mesh.vertices(:, 1), mesh.vertices(:, 2), mesh.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.5);
hold on;
quiver3(mesh.vertices(:, 1), mesh.vertices(:, 2), mesh.vertices(:, 3), mesh.normals(:, 1), mesh.normals(:, 2), mesh.normals(:, 3), 'Color', [0.8, 0.8, 0.8]);
hold off;
view([-90, 0]);
camlight;
view([90, 0]);
camlight;
axis equal;
grid off;
lighting gouraud;
axis off;
pause;

% Load spheres
file = fopen('mesh/new_centers.csv', 'r');
centers = textscan(file, '%f%f%f%[^\n\r]', 'Delimiter', ',', 'ReturnOnError', false);
centers = [centers{1 : end - 1}];
fclose(file);
load('mesh/radii.mat');
load('mesh/blocks.mat');
mesh.spheres = cell(1, size(centers, 1));
for i = 1 : length(mesh.spheres)
    % TODO rescale spheres
    mesh.spheres{i}.center = centers(i, :);
    mesh.spheres{i}.radius = radii{i};
end
mesh.blocks = blocks;
mesh.spheres{1}.bone = 18;
mesh.spheres{2}.bone = 17;
mesh.spheres{3}.bone = 16;
mesh.spheres{4}.bone = 15;
mesh.spheres{5}.bone = 14;
mesh.spheres{6}.bone = 13;
mesh.spheres{7}.bone = 12;
mesh.spheres{8}.bone = 2;
mesh.spheres{9}.bone = 11;
mesh.spheres{10}.bone = 10;
mesh.spheres{11}.bone = 9;
mesh.spheres{12}.bone = 2;
mesh.spheres{13}.bone = 8;
mesh.spheres{14}.bone = 7;
mesh.spheres{15}.bone = 6;
mesh.spheres{16}.bone = 2;
mesh.spheres{17}.bone = 5;
mesh.spheres{18}.bone = 4;
mesh.spheres{19}.bone = 3;
mesh.spheres{20}.bone = 2;
mesh.spheres{21}.bone = 2;
mesh.spheres{22}.bone = 1;
mesh.spheres{23}.bone = 1;
mesh.spheres{24}.bone = 1;
mesh.spheres{25}.bone = 1;
mesh.spheres{26}.bone = 1;
mesh.spheres{27}.bone = 1;
mesh.spheres{28}.bone = 2;
mesh.spheres{29}.bone = 2;
mesh.spheres{30}.bone = 2;

% Show spheres
v = get(gca, 'view');
trimesh(mesh.faces, mesh.vertices(:, 1), mesh.vertices(:, 2), mesh.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.5);
hold on;
for i = 1 : length(mesh.spheres)
    [X, Y, Z] = sphere;
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

% Load bones
mesh.bones = cell(1, 18);
mesh.bones{1}.parent = 0;
mesh.bones{2}.parent = 1;
mesh.bones{3}.parent = 2;
mesh.bones{4}.parent = 3;
mesh.bones{5}.parent = 4;
mesh.bones{6}.parent = 2;
mesh.bones{7}.parent = 6;
mesh.bones{8}.parent = 7;
mesh.bones{9}.parent = 2;
mesh.bones{10}.parent = 9;
mesh.bones{11}.parent = 10;
mesh.bones{12}.parent = 2;
mesh.bones{13}.parent = 12;
mesh.bones{14}.parent = 13;
mesh.bones{15}.parent = 2;
mesh.bones{16}.parent = 15;
mesh.bones{17}.parent = 16;
mesh.bones{18}.parent = 17;
mesh.bones{1}.center = mesh.spheres{26}.center * 0.5 + mesh.spheres{27}.center * 0.5;
mesh.bones{2}.center = mesh.spheres{22}.center;
mesh.bones{3}.center = mesh.spheres{20}.center;
mesh.bones{4}.center = mesh.spheres{19}.center;
mesh.bones{5}.center = mesh.spheres{18}.center;
mesh.bones{6}.center = mesh.spheres{16}.center;
mesh.bones{7}.center = mesh.spheres{15}.center;
mesh.bones{8}.center = mesh.spheres{14}.center;
mesh.bones{9}.center = mesh.spheres{12}.center;
mesh.bones{10}.center = mesh.spheres{11}.center;
mesh.bones{11}.center = mesh.spheres{10}.center;
mesh.bones{12}.center = mesh.spheres{8}.center;
mesh.bones{13}.center = mesh.spheres{7}.center;
mesh.bones{14}.center = mesh.spheres{6}.center;
mesh.bones{15}.center = mesh.spheres{8}.center * 0.5 + mesh.spheres{23}.center * 0.5;
mesh.bones{16}.center = mesh.spheres{4}.center;
mesh.bones{17}.center = mesh.spheres{3}.center;
mesh.bones{18}.center = mesh.spheres{2}.center;
axes = bone_axes(mesh.spheres);

% Show bone hierarchy
v = get(gca, 'view');
colors = hsv(length(mesh.bones));
trimesh(mesh.faces, mesh.vertices(:, 1), mesh.vertices(:, 2), mesh.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.5);
hold on;
for i = 1 : length(mesh.bones)
    [X, Y, Z] = sphere;
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

% Load weights
file = fopen('mesh/weights5.csv', 'r');
data = textscan(file, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]', 'Delimiter', ',', 'EmptyValue', NaN, 'HeaderLines', 1, 'ReturnOnError', false);
fclose(file);
mesh.weights = [data{1 : end - 1}];
mesh.weights = mesh.weights ./ repmat(sum(mesh.weights, 2), 1, size(mesh.weights, 2));

% Show weights
for i = 1 : length(mesh.bones)
    v = get(gca, 'view');
    trisurf(mesh.faces, mesh.vertices(:, 1), mesh.vertices(:, 2), mesh.vertices(:, 3), mesh.weights(:, i), 'EdgeColor', 'none');
    view([-90, 0]);
    camlight;
    view([90, 0]);
    camlight;
    axis equal;
    grid off;
    lighting gouraud;
    axis off;
    set(gca, 'view', v);
    colormap((0 : 0.01 : 1)' * colors(i, :));
    pause;
end

% Compute partitioning
[~, mesh.assignments] = max(mesh.weights, [], 2);
for i = 1 : length(mesh.bones)
    keep = mesh.assignments == i;
    indices = 1 : size(mesh.vertices, 1);
    mesh.bones{i}.indices = indices(keep);
    [mesh.bones{i}.vertices, mesh.bones{i}.faces] = mesh_filter_vertices(mesh.vertices, mesh.faces, keep);
    mesh.bones{i}.normals = mesh.normals(keep, :);
end

% Precompute HRBF
for i = 1 : length(mesh.bones)
    [mesh.bones{i}.hrbf.points, mesh.bones{i}.hrbf.normals] = hrbf_points(mesh.bones{i}.vertices, mesh.bones{i}.faces);
    mesh.bones{i}.hrbf.centers = hrbf_centers(mesh.bones{i}.hrbf.points);
    mesh.bones{i}.hrbf.coefficients = hrbf_coefficients(mesh.bones{i}.hrbf.centers, mesh.bones{i}.hrbf.points, mesh.bones{i}.hrbf.normals);
end

% Show HRBF reconstruction
for i = 1 : length(mesh.bones)
    trimesh(mesh.bones{i}.faces, mesh.bones{i}.vertices(:, 1), mesh.bones{i}.vertices(:, 2), mesh.bones{i}.vertices(:, 3), 'FaceAlpha', 0, 'EdgeColor', colors(i, :));
    view([-90, 0]);
    camlight;
    view([90, 0]);
    camlight;
    axis equal;
    grid off;
    lighting gouraud;
    axis off;
    hold on;
    [X, Y, Z] = meshgrid( ...
        min(mesh.bones{i}.vertices(:, 1)) - 2 : 0.4 : max(mesh.bones{i}.vertices(:, 1)) + 2, ...
        min(mesh.bones{i}.vertices(:, 2)) - 2 : 0.4 : max(mesh.bones{i}.vertices(:, 2)) + 2, ...
        min(mesh.bones{i}.vertices(:, 3)) - 2 : 0.4 : max(mesh.bones{i}.vertices(:, 3)) + 2 ...
    );
    V = hrbf_apply(mesh.bones{i}.hrbf.centers, mesh.bones{i}.hrbf.coefficients, [X(:), Y(:), Z(:)]);
    V = reshape(V, size(X));
    p = patch(isosurface(X, Y, Z, V, 0));
    isonormals(X, Y, Z, V, p);
    p.FaceColor = colors(i, :);
    p.FaceAlpha = 0.5;
    p.EdgeColor = 'none';
    quiver3(mesh.bones{i}.hrbf.points(:, 1), mesh.bones{i}.hrbf.points(:, 2), mesh.bones{i}.hrbf.points(:, 3), mesh.bones{i}.hrbf.normals(:, 1), mesh.bones{i}.hrbf.normals(:, 2), mesh.bones{i}.hrbf.normals(:, 3), 'Color', colors(i, :));
    hold off;
    pause;
end

% Show HRBF reconstruction (everything in same plot)
trimesh(mesh.faces, mesh.vertices(:, 1), mesh.vertices(:, 2), mesh.vertices(:, 3), 'FaceAlpha', 0, 'EdgeColor', [0.8, 0.8, 0.8]);
hold on;
for i = 1 : length(mesh.bones)
    [X, Y, Z] = meshgrid( ...
        min(mesh.bones{i}.vertices(:, 1)) - 2 : 0.4 : max(mesh.bones{i}.vertices(:, 1)) + 2, ...
        min(mesh.bones{i}.vertices(:, 2)) - 2 : 0.4 : max(mesh.bones{i}.vertices(:, 2)) + 2, ...
        min(mesh.bones{i}.vertices(:, 3)) - 2 : 0.4 : max(mesh.bones{i}.vertices(:, 3)) + 2 ...
    );
    V = hrbf_apply(mesh.bones{i}.hrbf.centers, mesh.bones{i}.hrbf.coefficients, [X(:), Y(:), Z(:)]);
    V = reshape(V, size(X));
    p = patch(isosurface(X, Y, Z, V, 0));
    isonormals(X, Y, Z, V, p);
    p.FaceColor = colors(i, :);
    p.FaceAlpha = 0.5;
    p.EdgeColor = 'none';
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
pause;

% Compute sphere SDF
% TODO update to last sphere model?
mesh.bones{1}.blocks = mesh.blocks([25, 26]);
mesh.bones{2}.blocks = mesh.blocks([18, 19, 20, 21, 22, 23, 24, 28, 29]);
mesh.bones{3}.blocks = mesh.blocks(15);
mesh.bones{4}.blocks = mesh.blocks(14);
mesh.bones{5}.blocks = mesh.blocks(13);
mesh.bones{6}.blocks = mesh.blocks(12);
mesh.bones{7}.blocks = mesh.blocks(11);
mesh.bones{8}.blocks = mesh.blocks(10);
mesh.bones{9}.blocks = mesh.blocks(9);
mesh.bones{10}.blocks = mesh.blocks(8);
mesh.bones{11}.blocks = mesh.blocks(7);
mesh.bones{12}.blocks = mesh.blocks(6);
mesh.bones{13}.blocks = mesh.blocks(5);
mesh.bones{14}.blocks = mesh.blocks(4);
mesh.bones{15}.blocks = mesh.blocks([16, 17, 27]);
mesh.bones{16}.blocks = mesh.blocks(3);
mesh.bones{17}.blocks = mesh.blocks(2);
mesh.bones{18}.blocks = mesh.blocks(1);

% Show sphere reconstruction
clf;
for i = 1 : length(mesh.bones)
    low = inf(1, 3);
    high = -inf(1, 3);
    for j = 1 : length(mesh.bones{i}.blocks)
        for k = mesh.bones{i}.blocks{j};
            s = mesh.spheres{k};
            low = min(low, s.center - s.radius - 1);
            high = max(high, s.center + s.radius + 1);
        end
    end
    [X, Y, Z] = meshgrid(low(1) : 0.8 : high(1), low(2) : 0.8 : high(2), low(3) : 0.8 : high(3));
    disp(['Computing sphere SDF for bone #', num2str(i), '...']);
    V = sphere_apply(mesh.spheres, mesh.bones{i}.blocks, [X(:), Y(:), Z(:)]);
    V = reshape(V, size(X));
    p = patch(isosurface(X, Y, Z, V, 0));
    isonormals(X, Y, Z, V, p);
    p.FaceColor = colors(i, :);
    p.FaceAlpha = 0.5;
    p.EdgeColor = 'none';
    if i == 2
        hold on;
    end
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
pause;

% Precompute implicit skinning related values
mesh.implicit.adj_mat = adjacency_matrix(mesh.faces);
[mesh.implicit.mvc, mesh.implicit.is_on_side] = mesh_mvc(mesh);

% Define union for palm, the rest uses bulge
blend = ones(1, length(mesh.bones));
for i = 1:length(mesh.bones)
    mesh.bones{i}.blend = 'bulge';
end
for i = [1, 2, 15]
    mesh.bones{i}.blend = 'union';
    blend(i) = 0;
end

mesh.implicit.blend_weight = zeros(size(mesh.vertices, 1), 1);
for i = 1 : size(mesh.vertices, 1)
    w = zeros(length(mesh.bones), 1);
    for b = 1 : length(mesh.bones)
        local = inv(axes{b}) * [mesh.vertices(i, :)'; 1];
        angle = atan2(local(3), local(2));
        w(b) = max(-sin(angle), 0);
    end
    mesh.implicit.blend_weight(i) = mesh.weights(i, :) * (w .* blend');
end
v = get(gca, 'view');
patch('Faces', mesh.faces, 'Vertices', mesh.vertices, 'FaceVertexCData', mesh.implicit.blend_weight * [1, 0.5, 0], ...
    'EdgeColor', 'none', 'FaceColor', 'interp', 'CDataMapping', 'scaled')
view([-90, 0]);
camlight;
view([90, 0]);
camlight;
axis equal;
grid off;
lighting gouraud;
axis off;
set(gca, 'view', v);
colormap hot;
pause;

% Save mesh
save('mesh/neutral.mat', 'mesh');
disp('neutral.mat saved');
