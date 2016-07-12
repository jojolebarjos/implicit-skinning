function [mean_value_coord, is_on_side] = mesh_mvc(mesh)
%MESH_MVC Compute the Mean value coordinates for each vertex with respect to its neighbours
%
% Inputs
%  mesh: complete mesh
%
% Outputs
%  mean_value_coord: n by n sparse double, The coefficients for the neighbours of v_i are in the column i
%

    [neighbours, is_on_side] = mesh_neighbours_list(mesh);

    points_count = size(mesh.vertices, 1);
    mean_value_coord = spalloc(points_count, points_count, 6 * points_count);
    for i = 1:points_count
        if ~is_on_side(i)
            grad = mesh.normals(i, :);
            % Projects the one-ring neighbours of v_i onto tangent plane
            ni = mesh.vertices(neighbours{i}, :);
            n_count = size(ni, 1);
            % q_proj = q - dot(q - p, n) * n, where n is a unit vector
            ni_proj = ni - dot(ni - repmat(mesh.vertices(i, :), [n_count, 1]), repmat(grad, [n_count, 1]), 2) * grad;

            mean_value_coord(neighbours{i}, i) = mvc(mesh.vertices(i, :), ni_proj)';
        end
    end
end

