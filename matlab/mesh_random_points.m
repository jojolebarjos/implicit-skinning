function [points, points_normals] = mesh_random_points(vertices, faces, count, poisson)
%MESH_RANDOM_POINTS Sample random points on mesh using Poisson sampling
%
% Inputs
%  vertices: v by 3 double
%  faces: f by 3 integer
%  count: integer (default v)
%  poisson: use Poisson sampling (default true)
%
% Outputs
%  points: count by 3 double
%  points_normals: count by 3 double
%

    % Define missing arguments
    if nargin < 3
        count = size(vertices, 1);
    end
    if nargin < 4
        poisson = true;
    end
    
    % Poisson sampling doesn't work if some vertices are unused
    [vertices, indices] = remove_unreferenced(vertices, faces);
    faces = indices(faces);
    
    % Generate points and normals
    if poisson
        [points, indices] = random_points_on_mesh(vertices, faces, count, 'Color', 'blue');
    else
        [points, indices] = random_points_on_mesh(vertices, faces, count, 'Color', 'white');
    end
    if nargout > 1
        points_normals = normals(vertices, faces);
        points_normals = normalizerow(points_normals(indices, :));
    end

end

