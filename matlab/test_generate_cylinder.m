
% Register gptoolbox
addpath(genpath('external'));

m = 100; n = 100;
[V,F] = create_regular_grid(m,n,1,0);
V = [sin(2*pi*V(:,1)), cos(2*pi*V(:,1)), (n-1)*10/(m-1)*V(:,2)];
mesh = {};
mesh.vertices = V;
mesh.faces = F;
mesh.normals = per_vertex_normals(V,F);
mesh.spheres{1}.center = [0, 0, 1]; mesh.spheres{1}.radius = 1; mesh.spheres{1}.bone = 1;
mesh.spheres{2}.center = [0, 0, 5]; mesh.spheres{2}.radius = 1; mesh.spheres{2}.bone = 1;
mesh.spheres{3}.center = [0, 0, 9]; mesh.spheres{3}.radius = 1; mesh.spheres{3}.bone = 2;
mesh.blocks{1} = [1, 2];
mesh.blocks{2} = [2, 3];
mesh.bones{1}.parent = 0; mesh.bones{1}.center = [0, 0, 1]; mesh.bones{1}.blocks{1} = [1, 2]; mesh.bones{1}.blend = 'bulge';
mesh.bones{2}.parent = 1; mesh.bones{2}.center = [0, 0, 5]; mesh.bones{2}.blocks{1} = [2, 3]; mesh.bones{2}.blend = 'bulge';
W = bone_heat(V, F, [0, 0, 1; 0, 0, 5; 0, 0, 9], 1:3, [1, 2; 2, 3], []);
mesh.weights = [W(:, 1)+W(:, 4), W(:, 3)+W(:, 5)];
mesh.implicit.adj_mat = adjacency_matrix(mesh.faces);
% TODO: Add HRBFs
[mesh.implicit.mvc, mesh.implicit.is_on_side] = mesh_mvc(mesh);
mesh.implicit.blend_weight = max(mesh.vertices(:, 2), zeros(size(mesh.vertices, 1), 1));

% Save mesh
save('mesh/cylinder.mat', 'mesh');
disp('cylinder.mat saved');
