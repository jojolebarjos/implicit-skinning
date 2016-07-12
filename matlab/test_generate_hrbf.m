
function test_generate_hrbf()

% Register gptoolbox
addpath(genpath('external'));

% Load neutral mesh
mesh = {};
load('mesh/neutral.mat');

% Regenerate and show HRBF
colors = hsv(length(mesh.bones));
i = 1;
offset = 1;
function regen(compute)
    
    % Recompute
    if compute
        [mesh.bones{i}.hrbf.points, mesh.bones{i}.hrbf.normals] = hrbf_points(mesh.bones{i}.vertices, mesh.bones{i}.faces, offset);
        mesh.bones{i}.hrbf.centers = hrbf_centers(mesh.bones{i}.hrbf.points);
        mesh.bones{i}.hrbf.coefficients = hrbf_coefficients(mesh.bones{i}.hrbf.centers, mesh.bones{i}.hrbf.points, mesh.bones{i}.hrbf.normals);
    end
    
    % Show
    v = get(gca, 'view');
    trimesh(mesh.bones{i}.faces, mesh.bones{i}.vertices(:, 1), mesh.bones{i}.vertices(:, 2), mesh.bones{i}.vertices(:, 3), 'FaceAlpha', 0, 'EdgeColor', colors(i, :));
    view([-90, 0]);
    camlight;
    view([90, 0]);
    camlight;
    axis equal;
    grid off;
    lighting gouraud;
    axis off;
    hold on;
    [X, Y, Z] = meshgrid( ...
        min(mesh.bones{i}.vertices(:, 1)) - 3 : 0.2 : max(mesh.bones{i}.vertices(:, 1)) + 3, ...
        min(mesh.bones{i}.vertices(:, 2)) - 3 : 0.2 : max(mesh.bones{i}.vertices(:, 2)) + 3, ...
        min(mesh.bones{i}.vertices(:, 3)) - 3 : 0.2 : max(mesh.bones{i}.vertices(:, 3)) + 3 ...
    );
    V = hrbf_apply(mesh.bones{i}.hrbf.centers, mesh.bones{i}.hrbf.coefficients, [X(:), Y(:), Z(:)]);
    V = reshape(V, size(X));
    p = patch(isosurface(X, Y, Z, V, 0));
    isonormals(X, Y, Z, V, p);
    p.FaceColor = colors(i, :);
    p.FaceAlpha = 0.5;
    p.EdgeColor = 'none';
    quiver3(mesh.bones{i}.hrbf.points(:, 1), mesh.bones{i}.hrbf.points(:, 2), mesh.bones{i}.hrbf.points(:, 3), mesh.bones{i}.hrbf.normals(:, 1), mesh.bones{i}.hrbf.normals(:, 2), mesh.bones{i}.hrbf.normals(:, 3), 'Color', colors(i, :));
    hold off;
    set(gca, 'view', v);
    
    % Handle key callback
    function callback(~, event, ~)
        switch event.Key
            case 's'
                i = mod(i, length(mesh.bones)) + 1;
                offset = 1;
                regen(0);
            case 'y'
                save('mesh/neutral.mat', 'mesh');
                disp('neutral.mat saved');
                i = mod(i, length(mesh.bones)) + 1;
                offset = 1;
                regen(0);
            case 'n'
                regen(1);
            case 'p'
                offset = offset + 0.25;
                regen(1);
            case 'm'
                offset = offset - 0.25;
                regen(1);
        end
    end
    set(gcf, 'KeyPressFcn', @callback);
end
regen(0);
pause;

end