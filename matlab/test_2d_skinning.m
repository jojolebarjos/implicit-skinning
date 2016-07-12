
% Load geometry
load('mesh/finger.mat');
c1 = [finger.bones{1}.center', 0];
r1 = finger.bones{1}.radius;
w1 = finger.W(:, 1);
c2 = [finger.bones{2}.center', 0];
r2 = finger.bones{2}.radius;
w2 = finger.W(:, 2);
c3 = [finger.bones{3}.center', 0];
r3 = finger.bones{3}.radius;
w3 = finger.W(:, 3);
v = [finger.x, finger.y, zeros(size(finger.x))];
N = length(finger.x);

% Configure mapping to have pure union at rest and maximal blend at 45°
map = @(alpha) sdf_blend_theta_curve(alpha, 0.1, pi / 4, pi / 2, pi / 4, 0, pi / 32, 1);
plot(linspace(0, pi), map(linspace(0, pi)));
pause;

% Approximate mesh using two bispheres
f1 = @(p) sdf_reparam(sphere_double(c1, r1, c2, r2, p), 1);
f2 = @(p) sdf_reparam(sphere_double(c2, r2, c3, r3, p), 1);
f = @(p) sdf_blend(f1(p), f2(p), map(0));

% Compute field and vertex offsets at rest
off = f(v);

% Compute barycentric coordinates (in 2D, this is easy, only 2 neighbours)
bar = zeros(N, 1);
for i = 2 : N-1
    p = v(i, :) - v(i + 1, :);
    d = v(i - 1, :) - v(i + 1, :);
    bar(i) = sum(p .* d) / (d * d');
end

% Simulate skinning for various angles
for alpha = [pi / 8, pi / 4, pi / 2, 3 * pi / 4]
    theta = map(alpha);
    
    % Reset geometry
    c1 = [finger.bones{1}.center', 0];
    c2 = [finger.bones{2}.center', 0];
    c3 = [finger.bones{3}.center' 0];
    v = [finger.x, finger.y, zeros(size(finger.x))];
    
    % Compute transform matrices
    m1 = eye(4);
    m2 = eye(4) * m1;
    m3 = matrix_rotation(-alpha, [0, 0, 1], c2) * m2;
    
    % Compute bone locations
    c1 = matrix_apply(m1, c1);
    c2 = matrix_apply(m2, c2);
    c3 = matrix_apply(m3, c3);
    
    % Compute initial guess using linear skinning (should use quaternion, though)
    for i = 1 : size(v, 1);
        m = m1 * w1(i) + m2 * w2(i) + m3 * w3(i);
        v(i, :) = matrix_apply(m, v(i, :));
    end
    
    % Recompute field
    f1 = @(p) sdf_reparam(sphere_double(c1, r1, c2, r2, p), 1);
    f2 = @(p) sdf_reparam(sphere_double(c2, r2, c3, r3, p), 1);
    f = @(p) sdf_blend(f1(p), f2(p), map(0));
    
    % Show result
    [X, Y] = meshgrid(linspace(-5, 5), linspace(-5, 5));
    Z = zeros(size(X));
    Z = reshape(f([X(:), Y(:), Z(:)]), size(X));
    contourf(X, Y, Z, 12);
    hold on;
    plot(v(:, 1), v(:, 2), 'r');
    hold off;
    axis square;
    title(['\alpha = ', num2str(alpha), ': initial guess']);
    pause;
    
    % Numerical derivative
    df = @(p) [ ...
        f([p(1) + 0.01, p(2), 0]) - f([p(1) - 0.01, p(2), 0]), ...
        f([p(1), p(2) + 0.01, 0]) - f([p(1), p(2) - 0.01, 0]) ...
    ] / 0.02;
    
    % Prepare self-intersection heuristic
    beta = zeros(N, 1);
    angle = zeros(N, 1);
    for i = 1 : N
        d = df(v(i, :));
        angle(i) = atan2(d(2), d(1));
    end

    % Move vertices
    for n = 1 : 10
        
        % Tangential relaxation
        if n > 1
            delta = zeros(N, 3);
            for i = 2 : N - 1
                if beta(i) == 0
                    p = v(i, :) - v(i + 1, :);
                    d = v(i - 1, :) - v(i + 1, :);
                    b = sum(p .* d) / (d * d');
                    mu = max(0, 1 - (abs(f(v(i, :)) - off(i)) - 1) ^ 4);
                    delta(i, :) = mu * d * (bar(i) - b);
                end
            end
            v = v + delta;
            
            % Show result
            [X, Y] = meshgrid(linspace(-5, 5), linspace(-5, 5));
            Z = zeros(size(X));
            Z = reshape(f([X(:), Y(:), Z(:)]), size(X));
            contourf(X, Y, Z, 12);
            hold on;
            plot(v(:, 1), v(:, 2), 'r');
            hold off;
            axis square;
            title(['\alpha = ', num2str(alpha), ': tangential relaxation #', num2str(n - 1)]);
            pause(0.1);
        end
        
        % Vertex projection
        sigma = 0.35;
        for i = 1 : N
            if beta(i) == 0
                d = [df(v(i, :)), 0];
                a = atan2(d(2), d(1));
                if mod(abs(a - angle(i)), 2 * pi) > 55 * pi / 180
                    beta(i) = 1;
                else
                    v(i, :) = v(i, :) - sigma * (f(v(i, :)) - off(i)) * d / sqrt(d * d');
                    angle(i) = a;
                end
            end
        end
        
        % Show result
        [X, Y] = meshgrid(linspace(-5, 5), linspace(-5, 5));
        Z = zeros(size(X));
        Z = reshape(f([X(:), Y(:), Z(:)]), size(X));
        contourf(X, Y, Z, 12);
        hold on;
        plot(v(:, 1), v(:, 2), 'r');
        hold off;
        axis square;
        title(['\alpha = ', num2str(alpha), ': vertex projection #', num2str(n)]);
        pause(0.1);
        
    end
    
    % Laplacian smoothing
    for i = 1 : 3
        beta = min(conv([0.1, 0.2, 1.0, 0.2, 0.1], beta), 1);
        beta = beta(3 : N+2);
    end
    centroid = [zeros(1, 3); [v(1 : N - 2, :); zeros(2, 3)] + [zeros(2, 3); v(3 : N, :)]; zeros(1, 3)] * 0.5;
    for i = 2 : N - 1
        v(i, :) = (1 - beta(i)) * v(i, :) + beta(i) * centroid(i, :);
    end
    
    % Show result
    [X, Y] = meshgrid(linspace(-5, 5), linspace(-5, 5));
    Z = zeros(size(X));
    Z = reshape(f([X(:), Y(:), Z(:)]), size(X));
    contourf(X, Y, Z, 12);
    hold on;
    plot(v(:, 1), v(:, 2), 'r');
    hold off;
    axis square;
    title(['\alpha = ', num2str(alpha), ': laplacian smoothing']);
    pause;
    
end
