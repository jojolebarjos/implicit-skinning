function centers = hrbf_centers(vertices, faces, count)
%HRBF_CENTERS Generate centers that should approximate given vertices
%
% Inputs
%  vertices: v by 3 double
%  faces: f by 3 integer
%  count: integer (default 10)
%
% Outputs
%  centers: count by 3 double
%

    % Define missing arguments
    if nargin < 2
        faces = convhull(vertices(:, 1), vertices(:, 2), vertices(:, 3));
    end
    if nargin < 3
        count = 10;
    end
    
    % Use Poisson sampling to select centers on mesh surface
    centers = mesh_random_points(vertices, faces, count);
    % TODO allow centers inside convex hull?

end

