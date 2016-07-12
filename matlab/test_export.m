
% Register gptoolbox
addpath(genpath('external'));

% Load mesh and poses
load('mesh/neutral.mat');
load('mesh/poses.mat');
%mesh = mesh_densify(mesh);

% Export mesh
file = fopen('mesh/neutral.obj', 'w');
for i = 1 : size(mesh.vertices, 1)
    fprintf(file, 'v %f %f %f\r\n', ...
        mesh.vertices(i, 1), mesh.vertices(i, 2), mesh.vertices(i, 3));
    fprintf(file, 'vn %f %f %f\r\n', ...
        mesh.normals(i, 1), mesh.normals(i, 2), mesh.normals(i, 3));
end
for i = 1 : size(mesh.faces, 1)
    fprintf(file, 'f %d/%d %d/%d %d/%d\r\n', ...
        mesh.faces(i, 1), mesh.faces(i, 1), mesh.faces(i, 2), mesh.faces(i, 2), mesh.faces(i, 3), mesh.faces(i, 3));
end
fclose(file);
disp('neutral.obj saved');

% Export skeleton
file = fopen('mesh/neutral.skl', 'w');
for i = 1 : 18
    c = mesh.bones{i}.center;
    fprintf(file, 'b %d %d %f %f %f\r\n', i, mesh.bones{i}.parent, c(1), c(2), c(3));
    for j = 1 : size(mesh.bones{i}.hrbf.centers, 1)
        p = mesh.bones{i}.hrbf.centers(j, :);
        c = mesh.bones{i}.hrbf.coefficients(j, :);
        fprintf(file, 'h %f %f %f %f %f %f %f\r\n', p(1), p(2), p(3), c(2), c(3), c(4), c(1));
    end
end
for i = 1 : length(poses)
    fprintf(file, 'p %d\r\n', i);
    transforms = bone_local_transforms(mesh, bone_transforms(mesh, poses{i}));
    for j = 1 : 18
        t = transforms{j};
        fprintf(file, 't %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\r\n', ...
            t(1, 1), t(1, 2), t(1, 3), t(1, 4), ...
            t(2, 1), t(2, 2), t(2, 3), t(2, 4), ...
            t(3, 1), t(3, 2), t(3, 3), t(3, 4), ...
            t(4, 1), t(4, 2), t(4, 3), t(4, 4));
    end
end
for i = 1 : size(mesh.weights, 1)
    [W, I] = sort(mesh.weights(i, :), 'descend');
    fprintf(file, 'w %d/%f %d/%f %d/%f %d/%f\r\n', ...
        I(1), W(1), I(2), W(2), I(3), W(3), I(4), W(4));
end
fclose(file);
disp('neutral.skl saved');
