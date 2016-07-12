function mesh = mesh_densify(mesh)
%MESH_DENSIFY Add a vertex a each face center
%
% Inputs
%  mesh: complete mesh (at least with vertices and faces)
%
% Outputs
%  mesh: complete mesh
%

    % Generate new vertices at each face center
    % TODO provide other interpolation scheme?
    [vertices, faces] = upsample(mesh.vertices, mesh.faces);
    old_vertices_indices = 1 : size(mesh.vertices, 1);
    new_vertices_indices = size(mesh.vertices, 1) + 1 : size(vertices, 1);
    
    adj_mat = adjacency_matrix(faces);
    % Retrieve the indices of the 2 old vertices surrounding every new vertex
    [row, ~] = find(adj_mat(old_vertices_indices, new_vertices_indices));
    indices = reshape(row, 2, size(row, 1) / 2)';
    
    mesh.vertices = vertices;
    mesh.faces = faces;
    
    % Interpolate normals
    if isfield(mesh, 'normals')
        normals = normalizerow( ...
            mesh.normals(indices(:, 1), :) + ...
            mesh.normals(indices(:, 2), :) ...
        );
        mesh.normals = [mesh.normals; normals];
    end

    % Interpolate weights
    if isfield(mesh, 'weights')
        weights = ( ...
            mesh.weights(indices(:, 1), :) + ...
            mesh.weights(indices(:, 2), :) ...
        );
        weights = weights ./ repmat(sum(weights, 2), 1, size(weights, 2));
        mesh.weights = [mesh.weights; weights];

        % Update assignments
        [~, mesh.assignments] = max(mesh.weights, [], 2);
    end
    
    % Update implicit precomputations
    if isfield(mesh, 'implicit')
        mesh.implicit.adj_mat = adj_mat;
        [mesh.implicit.mvc, mesh.implicit.is_on_side] = mesh_mvc(mesh);

        blend_weight = ( ...
            mesh.implicit.blend_weight(indices(:, 1), :) + ...
            mesh.implicit.blend_weight(indices(:, 2), :) ...
        ) ./ 2;
        mesh.implicit.blend_weight = [mesh.implicit.blend_weight; blend_weight];
    end
    
    % Update bones
    if isfield(mesh, 'bones')
        for i = 1 : length(mesh.bones)
            keep = mesh.assignments == i;
            indices = 1 : size(mesh.vertices, 1);
            mesh.bones{i}.indices = indices(keep);
            [mesh.bones{i}.vertices, mesh.bones{i}.faces] = mesh_filter_vertices(mesh.vertices, mesh.faces, keep);
            mesh.bones{i}.normals = mesh.normals(keep, :);
        end
    end
end

