function transforms = compute_transfer_transforms(mesh, centers, radii)
%COMPUTE_TRANSFER_TRANSFORMS Compute transform for each bone to fit given spheres

    % Compute local scale factors
    mesh.centers = zeros(27, 3);
    mesh.radii = zeros(1, 27);
    for i = 1 : 27
        mesh.centers(i, :) = mesh.spheres{i}.center;
        mesh.radii(i) = mesh.spheres{i}.radius;
    end
    dfactor = @(i, j) norm(mean(centers(j, :), 2) - mean(centers(i, :), 2)) / norm(mean(mesh.centers(j, :), 2) - mean(mesh.centers(i, :), 2));
    sfactor = @(i) mean(radii(i)) / mean(mesh.radii(i));
    % TODO maybe use smallest radius for scale (or average)
    scale = ones(18, 3);
    
    % Palm
    scale(1, 1) = dfactor([26, 27], 22);
    scale(1, 2) = 0.7 * dfactor(4, 16) + 0.3 * dfactor(20, 23); % TODO fix this
    scale(1, 3) = sfactor([24, 25, 26, 27]);
    scale(2, 1) = dfactor(22, 12);
    scale(2, 2 : 3) = dfactor(8, 16);
    scale(15, 1) = dfactor([8, 23], 4);
    scale(15, 2 : 3) = sfactor(8);
    
    % Thumb
    %scale(3, 1) = dfactor(20, 19); % TODO fix this
    scale(3, 1) = dfactor(19, 18);
    scale(3, 2 : 3) = sfactor(20);
    %scale(3, 2 : 3) = sfactor(19);
    scale(4, 1) = dfactor(19, 18);
    scale(4, 2 : 3) = sfactor(19);
    scale(5, 1) = dfactor(18, 17);
    scale(5, 2 : 3) = sfactor(18);
    
    % Index
    scale(6, 1) = dfactor(16, 15);
    scale(6, 2 : 3) = sfactor(16);
    scale(7, 1) = dfactor(15, 14);
    scale(7, 2 : 3) = sfactor(15);
    scale(8, 1) = dfactor(14, 13);
    scale(8, 2 : 3) = sfactor(14);
    
    % Middle
    scale(9, 1) = dfactor(12, 11);
    scale(9, 2 : 3) = sfactor(12);
    scale(10, 1) = dfactor(11, 10);
    scale(10, 2 : 3) = sfactor(11);
    scale(11, 1) = dfactor(10, 9);
    scale(11, 2 : 3) = sfactor(10);
    
    % Ring
    scale(12, 1) = dfactor(8, 7);
    scale(12, 2 : 3) = sfactor(8);
    scale(13, 1) = dfactor(7, 6);
    scale(13, 2 : 3) = sfactor(7);
    scale(14, 1) = dfactor(6, 5);
    scale(14, 2 : 3) = sfactor(6);
    
    % Pinky
    scale(16, 1) = dfactor(4, 3);
    scale(16, 2 : 3) = sfactor(4);
    scale(17, 1) = dfactor(3, 2);
    scale(17, 2 : 3) = sfactor(3);
    scale(18, 1) = dfactor(2, 1);
    scale(18, 2 : 3) = sfactor(2);
    
    % Compute global transforms
    A = compute_bone_axes(mesh.spheres);
    B = compute_bone_axes(centers);
    transforms = cell(1, 18);
    for i = 1 : 18
        transforms{i} = B{i} * diag([scale(i, :), 1]) * inv(A{i});
    end

end
