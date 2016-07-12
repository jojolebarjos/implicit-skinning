function [points, normals] = hrbf_points(vertices, faces, offset)
%HRBF_POINTS Generate control points for specified mesh
%
% Inputs
%  vertices: v by 3 double
%  faces: f by 3 integer
%
% Outputs
%  points: n by 3 double
%  normals: n by 3 double
%

    % TODO improve this using hole_fill
    
    % Add missing arguments
    if nargin < 3
        offset = 1;
    end

    % Sample mesh vertices
    [points, normals] = mesh_random_points(vertices, faces, floor(size(vertices, 1) / 10));
    
    % Add control points on holes
    holes = hole_find(faces);
    for i = 1 : length(holes)
        
        % Fill hole roughly
        [patch_vertices, patch_faces] = hole_fill(vertices, faces, holes{i}, false);
        keep = (1 : size(patch_vertices, 1)) > size(vertices, 1) + (size(patch_vertices, 1) - size(vertices, 1)) / 3 * 2;
        [patch_vertices, patch_faces] = mesh_filter_vertices(patch_vertices, patch_faces, keep);
        
        % Select random points and move them
        if size(patch_vertices, 1) < 11
            [patch_vertices, patch_faces] = upsample(patch_vertices, patch_faces);
        end
        [p, n] = mesh_random_points(patch_vertices, patch_faces, floor(size(patch_vertices, 1) / 5), false);
        p = p + offset * n;
        
        % Add them to control points
        points = [points; p];
        normals = [normals; n];
    end
    
end
