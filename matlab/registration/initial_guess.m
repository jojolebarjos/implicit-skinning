function transform = initial_guess(vertices)
%INITIAL_GUESS Normalize, center and rotate mesh vertices

    % Normalize data using mean and standard deviation
    center = mean(vertices);
    scale = max(std(vertices));
    transform = [eye(3) / scale, -center' / scale; 0, 0, 0, 1];
    
    % Try to find good orientation using PCA
    axes = princomp(apply_matrix(transform, vertices));
    transform = [axes', zeros(3, 1); 0, 0, 0, 1] * transform;

end

