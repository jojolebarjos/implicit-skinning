function [vertices, faces] = mesh_filter_vertices(vertices, faces, keep)
%MESH_FILTER_VERTICES Keep only specified vertices and discard unused elements
%
% Inputs
%  vertices: v by 3 double
%  faces: f by 3 integer
%  keep: v by 1 logical
%
% Outputs
%  vertices: v by 3 double
%  faces: f by 3 integer
%

    % Drop unused faces
    indices = 1 : size(vertices, 1);
    faces = faces(all(ismember(faces, indices(keep)), 2), :);
    
    % Remap vertices
    map = cumsum(keep);
    faces = map(faces);
    vertices = vertices(keep, :);

end

