
% http://ch.mathworks.com/help/matlab/ref/uicontrol.html
% http://ch.mathworks.com/matlabcentral/answers/60745-set-x-y-z-data-in-trimesh-without-re-plotting

function test_spheres_manual()
    
    % Load data
    mesh = load('../mesh/neutral.mat');
    mesh = mesh.mesh;
    points = load('anastasia_fixed.mat');
    
    % Precompute guesses
    mesh.centers = zeros(27, 3);
    mesh.radii = zeros(1, 27);
    for i = 1 : 27
        mesh.centers(i, :) = mesh.spheres{i}.center;
        mesh.radii(i) = mesh.spheres{i}.radius;
    end
    dfactor = @(i, j) norm(mean(points.centers(j, :), 2) - mean(points.centers(i, :), 2)) / norm(mean(mesh.centers(j, :), 2) - mean(mesh.centers(i, :), 2));
    scale = ones(18, 1);
    scale(1) = dfactor([26, 27], 22);
    scale(2) = dfactor(22, 12);
    scale(15) = dfactor([8, 23], 4);
    scale(3) = dfactor(19, 18);
    scale(4) = dfactor(19, 18);
    scale(5) = dfactor(18, 17);
    scale(6) = dfactor(16, 15);
    scale(7) = dfactor(15, 14);
    scale(8) = dfactor(14, 13);
    scale(9) = dfactor(12, 11);
    scale(10) = dfactor(11, 10);
    scale(11) = dfactor(10, 9);
    scale(12) = dfactor(8, 7);
    scale(13) = dfactor(7, 6);
    scale(14) = dfactor(6, 5);
    scale(16) = dfactor(4, 3);
    scale(17) = dfactor(3, 2);
    scale(18) = dfactor(2, 1);
    % TODO clamp excessive values
    scale = [scale scale scale];
    guess = mean(scale(:));
    scale = scale / guess;
    scale = clamp(scale, 0.25, 4);
    A = compute_bone_axes(mesh.spheres);
    B = compute_bone_axes(points.centers);
    transforms = cell(1, 18);
    for i = 1 : 18
        transforms{i} = B{i} * diag([guess * scale(i, :), 1]) * inv(A{i});
    end
    
    % Hide figure until ready
    f = figure('Visible', 'off');
    
    % Show mesh
    modified = skin_linear(mesh, transforms);
    trimesh(points.faces, points.vertices(:, 1), points.vertices(:, 2), points.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.5);
    hold on;
    h = trimesh(modified.faces, modified.vertices(:, 1), modified.vertices(:, 2), modified.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.9, 0.4], 'FaceAlpha', 0.5);
    hold off;
    view([-90, 0]);
    camlight;
    view([90, 0]);
    camlight;
    axis equal;
    grid off;
    lighting gouraud;
    axis off;
    
    % Add controls
    for i = 1 : 18
        for j = 1 : 3
            uicontrol('Style', 'slider', 'Min', 0.25, 'Max', 4, 'Value', scale(i, j), 'Position', [j * 60 - 40, 25 * (19 - i) + 25, 50, 20], 'Callback', @update, 'Tag', num2str((i - 1) * 3 + j));
        end
    end
    uicontrol('Style', 'pushbutton', 'String', 'Save', 'Position', [20 20 50 20], 'Callback', @save_scale);
    uicontrol('Style', 'pushbutton', 'String', 'Load', 'Position', [90 20 50 20], 'Callback', @load_scale);
    
    % Show figure
    set(f, 'Visible', 'on')
    
    % Sliders callback
    function update(source, callbackdata)
        
        % Update transform
        modifier = get(source, 'Value');
        tag = str2num(get(source, 'Tag'));
        axe = mod(tag - 1, 3) + 1;
        bone = (tag - axe) / 3 + 1;
        scale(bone, axe) = modifier;
        transforms{bone} = B{bone} * diag([guess * scale(bone, :), 1]) * inv(A{bone});
        
        % Update mesh
        M = skin_linear(mesh, transforms);
        hold on;
        delete(h);
        h = trimesh(M.faces, M.vertices(:, 1), M.vertices(:, 2), M.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.9, 0.4], 'FaceAlpha', 0.5);
        lighting gouraud;
        hold off;
        
    end

    % Save scale
    function save_scale(source, callbackdata)
        save('scale.mat', 'scale');
    end

    % Load scale
    function load_scale(source, callbackdata)
        scale = load('scale.mat');
        scale = scale.scale;
        for i = 1 : 18
            transforms{i} = B{i} * diag([guess * scale(i, :), 1]) * inv(A{i});
        end
        M = skin_linear(mesh, transforms);
        hold on;
        delete(h);
        h = trimesh(M.faces, M.vertices(:, 1), M.vertices(:, 2), M.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.9, 0.4], 'FaceAlpha', 0.5);
        lighting gouraud;
        hold off;
    end

end


