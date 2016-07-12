
% Register gptoolbox
addpath(genpath('external'));

% Load neutral mesh and keep only a finger
load('mesh/neutral.mat');
distances = heat_geodesic(mesh.vertices, mesh.faces, 1);
distances = max(distances) - distances;
[finger.vertices, finger.faces] = mesh_filter_vertices(mesh.vertices, mesh.faces, distances < 2);
finger.normals = per_vertex_normals(finger.vertices, finger.faces);

% For some parts, test HRBF
phalange = mesh.bones{13};
palm = mesh.bones{2};
proto = mesh.bones{15};
parts = {finger, phalange, palm, proto};
for n = 1 : length(parts)
    part = parts{n};

    % Show mesh
    trimesh(part.faces, part.vertices(:, 1), part.vertices(:, 2), part.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8]);
    hold on;
    quiver3(part.vertices(:, 1), part.vertices(:, 2), part.vertices(:, 3), part.normals(:, 1), part.normals(:, 2), part.normals(:, 3));
    hold off;
    view([0, 90]);
    camlight;
    view([180, -90]);
    camlight;
    axis equal;
    grid off;
    lighting gouraud;
    axis off;
    pause;

    % Select centers and points and compute coefficients
    centers = hrbf_centers(part.vertices, part.faces);
    [points, normals] = hrbf_points(part.vertices, part.faces);

    % Show HRBF control points
    trimesh(part.faces, part.vertices(:, 1), part.vertices(:, 2), part.vertices(:, 3), 'EdgeColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0);
    hold on;
    for i = 1 : size(centers, 1)
        [X, Y, Z] = sphere;
        X = X / 4 + centers(i, 1);
        Y = Y / 4 + centers(i, 2);
        Z = Z / 4 + centers(i, 3);
        surf(X, Y, Z, 'FaceColor', 'black');
    end
    quiver3(points(:, 1), points(:, 2), points(:, 3), normals(:, 1), normals(:, 2), normals(:, 3), 'Color', [0.8, 0.5, 0.5]);
    hold off;
    view([0, 90]);
    camlight;
    view([180, -90]);
    camlight;
    axis equal;
    grid off;
    lighting gouraud;
    axis off;
    pause;

    % Compute coefficients
    coefficients = hrbf_coefficients(centers, points, normals);

    % Show reconstruction
    trimesh(part.faces, part.vertices(:, 1), part.vertices(:, 2), part.vertices(:, 3), 'EdgeColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0);
    hold on;
    [~, gradients] = hrbf_apply(centers, coefficients, part.vertices);
    quiver3(part.vertices(:, 1), part.vertices(:, 2), part.vertices(:, 3), gradients(:, 1), gradients(:, 2), gradients(:, 3), 'Color', [1.0, 0.5, 0.5]);
    [X, Y, Z] = meshgrid( ...
        min(part.vertices(:, 1)) - 2 : 0.2 : max(part.vertices(:, 1)) + 2, ...
        min(part.vertices(:, 2)) - 2 : 0.2 : max(part.vertices(:, 2)) + 2, ...
        min(part.vertices(:, 3)) - 2 : 0.2 : max(part.vertices(:, 3)) + 2 ...
    );
    V = hrbf_apply(centers, coefficients, [X(:), Y(:), Z(:)]);
    V = reshape(V, size(X));
    p = patch(isosurface(X, Y, Z, V, 0));
    isonormals(X, Y, Z, V, p);
    p.FaceColor = [0.8, 0.5, 0.5];
    p.FaceAlpha = 0.5;
    p.EdgeColor = 'none';
    for s = [-0.2, 0.4, 0.8]
        p = patch(isosurface(X, Y, Z, V, s));
        isonormals(X, Y, Z, V, p);
        p.FaceColor = [0.6, 0.5, 0.5];
        p.FaceAlpha = 0.1;
        p.EdgeColor = 'none';
    end
    hold off;
    view([0, 90]);
    camlight;
    view([180, -90]);
    camlight;
    axis equal;
    grid off;
    lighting gouraud;
    axis off;
    pause;

end
