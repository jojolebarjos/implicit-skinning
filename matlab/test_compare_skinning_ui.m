function test_compare_skinning_ui()

    % Register gptoolbox
    addpath(genpath('external'));

    % Load mesh
    mesh = load('mesh/neutral.mat');
    mesh = mesh.mesh;
    
    % Compute axis and initial transforms
    transforms = cell(1, 18);
    for i = 1 : 18
        transforms{i} = eye(4);
    end
    axes = bone_axes(mesh.spheres);

    % Prepare components
    f = figure('Visible', 'off');
    popup = uicontrol( ...
        'Style', 'popup', ...
        'String', {'linear', 'dualquat', 'implicit'}, ...
        'Position', [20 20 100 20], ...
        'Callback', @update ...
    );
    sliders = cell(1, 3);
    for i = 1 : 3
        sliders{i} = uicontrol( ...
            'Style', 'slider', ...
            'Min', -1, 'Max', 2, 'Value', 0, ...
            'Position', [20, 25 * (4 - i) + 20, 100, 20], ...
            'Callback', @update ...
        );
    end
    
    % Render mesh
    h = trimesh(mesh.faces, mesh.vertices(:, 1), mesh.vertices(:, 2), mesh.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.5);
    view([-90, 0]);
    camlight;
    view([90, 0]);
    camlight;
    axis equal;
    grid off;
    lighting gouraud;
    axis off;
    
    % Show figure
    set(f, 'Visible', 'on')
    
    % Update callback
    function update(varargin)
        transforms{6} = matrix_rotation( ...
            get(sliders{1}, 'Value'), ...
            matrix_apply(transforms{2}, axes{6}(1 : 3, 2)', 0), ...
            matrix_apply(transforms{2}, axes{6}(1 : 3, 4)') ...
        ) * transforms{2};
        transforms{7} = matrix_rotation( ...
            get(sliders{2}, 'Value'), ...
            matrix_apply(transforms{6}, axes{7}(1 : 3, 2)', 0), ...
            matrix_apply(transforms{6}, axes{7}(1 : 3, 4)') ...
        ) * transforms{6};
        transforms{8} = matrix_rotation( ...
            get(sliders{3}, 'Value'), ...
            matrix_apply(transforms{7}, axes{8}(1 : 3, 2)', 0), ...
            matrix_apply(transforms{7}, axes{8}(1 : 3, 4)') ...
        ) * transforms{7};
        transformed = mesh;
        switch get(popup, 'Value')
        case 1
            transformed = skin_linear(transformed, transforms);
        case 2
            transformed = skin_dualquat(transformed, transforms);
        case 3
            transformed = skin_implicit(transformed, transforms, 20, 'hrbf');
        end
        hold on;
        delete(h);
        h = trimesh(transformed.faces, transformed.vertices(:, 1), transformed.vertices(:, 2), transformed.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.5);
        lighting gouraud;
        hold off;
    end

end

