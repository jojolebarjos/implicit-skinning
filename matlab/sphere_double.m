function [d, g] = sphere_double(c1, r1, c2, r2, p)
%SPHERE_DOUBLE Compute distance and gradient of double sphere SDF
%
% Inputs
%  c1: 1 by 3 double
%  r1: double
%  c2: 1 by 3 double
%  r2: double
%  p: n by 3 double
%
% Outputs
%  d: n by 1 double
%  g: n by 3 double
%

    % First sphere cannot be smaller
    if r1 < r2
        [d, g] = sphere_double(c2, r2, c1, r1, p);
        return;
    end
    
    % Get main axe
    x = c2 - c1;
    x_n = norm(x);
    x = x / x_n;
    
    % Second sphere cannot completely be inside the first one
    if x_n <= r1 - r2
        [d, g] = sphere_single(c1, r1, p);
        return;
    end

    % Vectorize code
    n = size(p, 1);
    d = zeros(n, 1);
    g = zeros(n, 3);
    
    % Project point on axe
    c1p = p - repmat(c1, n, 1);
    q = (c1p * x') * x + repmat(c1, n, 1);
    
    % Get main normal
    y = p - q;
    y_n = sqrt(sum(y .* y, 2));
    y_n(abs(y_n) < 0.00001) = 1;
    y = y ./ repmat(y_n, 1, 3);
    
    % Compute slope
    e = (r1 - r2) ^ 2 / x_n;
    f = sqrt((x_n - e) * e);
    s = f / (x_n - e);
    
    % Compute tangential axes
    u = normalizerow(repmat(x, n, 1) - s * y);
    v = normalizerow(y + s * repmat(x, n, 1));
    
    % Project points on tangent
    alpha = dot(u, c1p, 2) / sqrt(x_n * x_n - (r1 - r2) ^ 2);
    
    % First sphere
    first = alpha <= 0;
    [d(first), g(first, :)] = sphere_single(c1, r1, p(first, :));
    
    % Second sphere
    second = alpha >= 1;
    [d(second), g(second, :)] = sphere_single(c2, r2, p(second, :));
    
    % In between
    between = ~second & ~first;
    d(between) = dot(c1p(between, :), v(between, :), 2) - r1;
    g(between, :) = v(between, :);

end
