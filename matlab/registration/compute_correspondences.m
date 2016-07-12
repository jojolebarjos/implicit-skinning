function pairs = compute_correspondences(mesh_vertices, mesh_normals, points_vertices, points_normals, distance_threshold, cos_angle_threshold)
%COMPUTE_CORRESPONDENCES For each mesh vertex, find closest point in cloud

    % Define missing arguments
    if nargin < 5
        distance_threshold = 0.1;
    end
    if nargin < 6
        cos_angle_threshold = cos(30 * pi / 180);
    end

    % For each mesh vertex, find closest point
    % TODO this is naive... use Kd-tree?
    n = size(mesh_vertices, 1);
    m = size(points_vertices, 1);
    pairs = zeros(n, 2);
    for i = 1 : n
        delta = points_vertices - repmat(mesh_vertices(i, :), m, 1);
        distances = sum(delta .^ 2, 2);
        [~, j] = min(distances);
        pairs(i, 1) = i;
        pairs(i, 2) = j;
    end
    
    % Filter bad pairs
    mesh_normals = normalizerow(mesh_normals);
    points_normals = normalizerow(points_normals);    
    % TODO vectorize this
    keep = true(n, 1);
    for i = 1 : n
        p = pairs(i, 1);
        q = pairs(i, 2);
        
        % Check if distance is too large
        distance = norm(mesh_vertices(p, :) - points_vertices(q, :));
        if distance > distance_threshold
            keep(i) = false;
        end
        
        % Check if normals differ
        cos_angle = mesh_normals(p, :) * points_normals(q, :)';
        if cos_angle < cos_angle_threshold
            keep(i) = false;
        end
        
    end
    pairs = pairs(keep, :);

end
