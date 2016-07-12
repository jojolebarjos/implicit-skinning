
% Register gptoolbox
addpath(genpath('external'));

% Load neutral mesh
load('mesh/neutral.mat');

% Fill holes of partitions
for i = 1 : length(mesh.bones)
    vertices = mesh.bones{i}.vertices;
    faces = mesh.bones{i}.faces;
    holes = hole_find(faces);
    count = size(vertices, 1);
    for j = 1 : length(holes)
        [vertices, faces] = hole_fill(vertices, faces, holes{j});
    end
    normals = per_vertex_normals(vertices, faces);
    colormap([0.8, 0.8, 0.8; 0.8, 0.5, 0.5]);
    colors = [ones(count, 1); 2 * ones(size(vertices, 1) - count, 1)];
    trimesh(faces, vertices(:, 1), vertices(:, 2), vertices(:, 3), colors, 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.5);
    view([-90, 0]);
    camlight;
    view([90, 0]);
    camlight;
    axis equal;
    grid off;
    lighting gouraud;
    axis off;
    pause;
end
