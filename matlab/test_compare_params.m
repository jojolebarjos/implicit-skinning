function test_compare_params()

    % Register gptoolbox
    addpath(genpath('external'));

    % Load mesh
    mesh = load('mesh/neutral.mat');
    mesh = mesh_densify(mesh.mesh);
%     mesh.implicit.blend_weight(mesh.implicit.blend_weight > 0.01) = 1;
%     mesh.implicit.blend_weight = mesh.implicit.blend_weight .* 2;
%     mesh.implicit.blend_weight = ones(size(mesh.vertices, 1), 1);
    reparam_threshold = 2;
    
    transforms = cellfun(@double, repmat({eye(4)}, 1, length(mesh.bones)), 'Un', 0);
    transforms{7} = matrix_rotation(pi / 3, [0, 0, 1], mesh.bones{7}.center);
    transforms{8} = transforms{7} * matrix_rotation(pi / 2, [0, 0, 1], mesh.bones{8}.center);
    
    filter = mesh.weights(:, 6) > 0.5 | mesh.weights(:, 7) ~= 0 | mesh.weights(:, 8) ~= 0;

    % Prepare components
    f = figure('Visible', 'off');
    sliders = cell(1, 8);
    values = [pi / 8, pi / 2, 1.1 * pi / 2, pi / 4, pi / 10, pi / 4, 0.5, 0.5];
    for i = 1 : 3
        sliders{i} = uicontrol( ...
            'Style', 'slider', ...
            'Min', 0, 'Max', pi, 'Value', values(i), ...
            'Position', [20, 25 * (9 - i), 100, 20], ...
            'Callback', @update ...
        );
    end
    for i = 4 : 6
        sliders{i} = uicontrol( ...
            'Style', 'slider', ...
            'Min', 0, 'Max', pi / 4, 'Value', values(i), ...
            'Position', [20, 25 * (9 - i), 100, 20], ...
            'Callback', @update ...
        );
    end
    for i = 7 : 8
        sliders{i} = uicontrol( ...
            'Style', 'slider', ...
            'Min', 0.1, 'Max', 10, 'Value', values(i), ...
            'Position', [20, 25 * (9 - i), 100, 20], ...
            'Callback', @update ...
        );
    end
    
    % Render mesh
    implicit = skin_implicit(mesh, transforms, reparam_threshold, 'hrbf');
    [V, F] = mesh_filter_vertices(implicit.vertices, implicit.faces, filter);
%     h = trimesh(F, V(:, 1), V(:, 2), V(:, 3), 'EdgeColor', [0.4, 0.25, 0.0], 'EdgeAlpha', 0.2, 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 1.0);
    h = trimesh(F, V(:, 1), V(:, 2), V(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 1.0, 'SpecularStrength', 0.4);
    l1 = [0, 30]; l2 = [180, 30]; vp = [90, 90];
    view(l1); camlight;
    view(l2); camlight;
    lighting gouraud;
    view(vp);
    axis equal;
    grid off;
    axis off;
    
    % Show figure
    set(f, 'Visible', 'on')
    
    % Update callback
    function update(varargin)
        mesh.implicit.alpha_to_theta.alpha0 = get(sliders{1}, 'Value');
        mesh.implicit.alpha_to_theta.alpha1 = get(sliders{2}, 'Value');
        mesh.implicit.alpha_to_theta.alpha2 = get(sliders{3}, 'Value');
        mesh.implicit.alpha_to_theta.theta0 = get(sliders{4}, 'Value');
        mesh.implicit.alpha_to_theta.theta1 = get(sliders{5}, 'Value');
        mesh.implicit.alpha_to_theta.theta2 = get(sliders{6}, 'Value');
        mesh.implicit.alpha_to_theta.omega0 = get(sliders{7}, 'Value');
        mesh.implicit.alpha_to_theta.omega1 = get(sliders{8}, 'Value');
        disp([mesh.implicit.alpha_to_theta.alpha0, ...
                mesh.implicit.alpha_to_theta.alpha1, ...
                mesh.implicit.alpha_to_theta.alpha2, ...
                mesh.implicit.alpha_to_theta.theta0, ...
                mesh.implicit.alpha_to_theta.theta1, ...
                mesh.implicit.alpha_to_theta.theta2, ...
                mesh.implicit.alpha_to_theta.omega0, ...
                mesh.implicit.alpha_to_theta.omega1]);
        implicit = skin_implicit(mesh, transforms, reparam_threshold, 'hrbf');
        hold on;
        delete(h);
        [V, F] = mesh_filter_vertices(implicit.vertices, implicit.faces, filter);
%         h = trimesh(F, V(:, 1), V(:, 2), V(:, 3), 'EdgeColor', [0.4, 0.25, 0.0], 'EdgeAlpha', 0.2, 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 1.0);
        h = trimesh(F, V(:, 1), V(:, 2), V(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 1.0, 'SpecularStrength', 0.4);
        lighting gouraud;
        hold off;
    end

end

