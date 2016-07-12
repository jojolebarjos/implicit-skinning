
% Load mesh
c1 = [-2.5, 0.5, 0];
r1 = 0.8;
c2 = [0, 0.6, 0];
r2 = 0.6;
c3 = [1.8, -1, 0];
r3 = 0.5;
c1c2 = c2 - c1;
c2c3 = c3 - c2;
alpha = acos(sum(c1c2 .* c2c3) / (norm(c1c2) * norm(c2c3)));
n12 = cross([0, 0, 1], c2 - c1);
n23 = cross([0, 0, 1], c3 - c2);

% Configure mapping to have pure union at rest and maximal blend at 45°
map = @(alpha) sdf_blend_theta_curve(alpha, 0.1, pi / 4, pi / 2, pi / 4, 0, pi / 32, 1);
plot(linspace(0, pi), map(linspace(0, pi)));
pause;

% Reference line, using pure union
f1 = @(p) sdf_reparam(sphere_double(c1, r1, c2, r2, p), 3);
f2 = @(p) sdf_reparam(sphere_double(c2, r2, c3, r3, p), 3);
f = @(p) sdf_union(f1(p), f2(p));
xs = linspace(0, 1);
ys = linspace(0, 1);
[X, Y] = meshgrid(xs, ys);
Z = zeros(size(X));
Z = f([X(:) * 10 - 5.5, Y(:) * 10 - 5, Z(:)]);
Z = reshape(Z, size(X));
contourf(X, Y, Z, 16);
hold on;
C = contourc(xs, ys, Z, [0.5, 0.5]);
C = C(1 : 2, 2 : C(2, 1) + 1);
plot(C(1, :), C(2, :), 'r', 'LineWidth', 2);
hold off;
axis square;
pause;

% Constant theta, based on bone angle
theta = map(alpha);
for threshold = 1:5
    f1 = @(p) sdf_reparam(sphere_double(c1, r1, c2, r2, p), threshold);
    f2 = @(p) sdf_reparam(sphere_double(c2, r2, c3, r3, p), threshold);
    f = @(p) sdf_blend(f1(p), f2(p), theta);
    xs = linspace(0, 1);
    ys = linspace(0, 1);
    [X, Y] = meshgrid(xs, ys);
    Z = zeros(size(X));
    Z = f([X(:) * 10 - 5.5, Y(:) * 10 - 5, Z(:)]);
    Z = reshape(Z, size(X));
    contourf(X, Y, Z, 16);
    hold on;
    C = contourc(xs, ys, Z, [0.5, 0.5]);
    C = C(1 : 2, 2 : C(2, 1) + 1);
    plot(C(1, :), C(2, :), 'r', 'LineWidth', 2);
    hold off;
    axis square;
    pause;
end

% Constant theta, based on bone angle, using blend weights
theta = map(alpha);
for threshold = 1:5
    w = @(p) 1 - ((p - repmat(c1, size(p, 1), 1)) * n12' > 0 | (p - repmat(c2, size(p, 1), 1)) * n23' > 0);
    f1 = @(p) sdf_reparam(sphere_double(c1, r1, c2, r2, p), threshold);
    f2 = @(p) sdf_reparam(sphere_double(c2, r2, c3, r3, p), threshold);
    f = @(p) w(p) .* sdf_blend(f1(p), f2(p), theta) + (1 - w(p)) .* sdf_union(f1(p), f2(p));
    xs = linspace(0, 1);
    ys = linspace(0, 1);
    [X, Y] = meshgrid(xs, ys);
    Z = zeros(size(X));
    Z = f([X(:) * 10 - 5.5, Y(:) * 10 - 5, Z(:)]);
    Z = reshape(Z, size(X));
    contourf(X, Y, Z, 16);
    hold on;
    C = contourc(xs, ys, Z, [0.5, 0.5]);
    C = C(1 : 2, 2 : C(2, 1) + 1);
    plot(C(1, :), C(2, :), 'r', 'LineWidth', 2);
    hold off;
    axis square;
    pause;
end

% Theta depends on angle between gradients
for threshold = 1:5
    xs = linspace(0, 1);
    ys = linspace(0, 1);
    [X, Y] = meshgrid(xs, ys);
    Z = zeros(size(X));
    P = [X(:) * 10 - 5.5, Y(:) * 10 - 5, Z(:)];
    [F1, G1] = sphere_double(c1, r1, c2, r2, P);
    [F2, G2] = sphere_double(c2, r2, c3, r3, P);
    F1 = sdf_reparam(F1, threshold);
    F2 = sdf_reparam(F2, threshold);
    G1 = -normalizerow(G1);
    G2 = -normalizerow(G2);
    A = acos(dot(G1, G2, 2));
    T = map(A);
    Z = sdf_blend(F1, F2, T);
    Z = reshape(Z, size(X));
    contourf(X, Y, Z, 16);
    hold on;
    C = contourc(xs, ys, Z, [0.5, 0.5]);
    C = C(1 : 2, 2 : C(2, 1) + 1);
    plot(C(1, :), C(2, :), 'r', 'LineWidth', 2);
    hold off;
    axis square;
    pause;
end
