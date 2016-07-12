function [neighbours, is_vertex_on_side] = mesh_neighbours_list(mesh)
%MESH_NEIGHBOURS_LIST For every vertex, compute the list of (one-ring) neighbours in counter-clockwise order
%
% Inputs
%  F: f by 3 double
%
% Outputs
%  neighbours: n by 1 array of x by 1 double, where x is the number of neighbours for the vertex
%  is_vertex_on_side: n by 1 logical, Specify if a vertex is on a mesh boundary
%

    is_on_side_indices = outline(mesh.faces);
    is_vertex_on_side = zeros(size(mesh.vertices, 1), 1);
    is_vertex_on_side(is_on_side_indices(:, 1)) = 1;

    neighbours = cell(size(mesh.vertices, 1), 1);
    for i = 1:size(mesh.vertices, 1)
        if ~is_vertex_on_side(i)
            v = mesh.vertices(i, :);
            neighbours_indices = find(mesh.implicit.adj_mat(:, i));
            neighbours_positions = mesh.vertices(neighbours_indices, :);

            % Compute unit vectors from v to the neighbours
            v_to_neighbours = neighbours_positions - repmat(v, [size(neighbours_positions, 1), 1]);
            v_to_neighbours = normalizerow(v_to_neighbours);
            
            % Get a 2d basis
            x_axis = v_to_neighbours(1, :);
            y_axis = cross(mesh.normals(i, :), x_axis);
            y_axis = y_axis ./ norm(y_axis);
            
            % Project the neighbours onto this bases
            x_coords = dot(repmat(x_axis, [size(v_to_neighbours, 1), 1]), v_to_neighbours, 2);
            y_coords = dot(repmat(y_axis, [size(v_to_neighbours, 1), 1]), v_to_neighbours, 2);
            
            % Get the angles with respect to the first neighbours around the normal vector
            angles = atan2(y_coords, x_coords);
            
            angles_with_indices = [angles, neighbours_indices];
            angles_with_indices = sortrows(angles_with_indices);
            
            neighbours{i} = angles_with_indices(:, 2);
        end
    end
end