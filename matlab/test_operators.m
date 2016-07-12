
N = 100;
[f1, f2] = meshgrid( ...
    linspace(0, 1, N), ...
    linspace(0, 1, N) ...
);

% Union (max)
g = reshape(sdf_union(f1(:), f2(:)), size(f1));
contourf(f1, f2, g, 32);
axis square;
hold on;
plot([0, 1], [0, 1], 'r');
hold off;
title('Union');
pause;

% Implicit skinning blend
thetas = [pi / 16, pi / 8, pi / 8 + pi / 16, pi / 4 - 0.05, pi / 4]; % between 0 and pi/4
for theta = thetas
    g = reshape(sdf_blend(f1(:), f2(:), theta), size(f1));
    contourf(f1, f2, g, 32);
    axis square;
    hold on;
    k1 = @(C) tan(theta) .* C;
    k2 = @(C) 1 / 2 * (7 - 5 * tan(theta)) * C + 7 / 4 * (tan(theta) - 1);
    x = [0, 1];
    plot(x, k1(x), 'r');
    plot(k1(x), x, 'r');
    plot(x, k2(x), 'r');
    plot(k2(x), x, 'r');
    axis([0, 1, 0, 1]);
    hold off;
    title(['Blend (theta=', num2str(theta), ')']);
    pause;
end

% Clean union
g = reshape(sdf_clean_union(f1(:), f2(:)), size(f1));
contourf(f1, f2, g, 32);
axis square;
hold on;
x = linspace(0, 1);
plot(x, sqrt(x / 2), 'r');
x = linspace(0, 1 / sqrt(2));
plot(x, 2 * x .* x, 'r');
hold off;
title('Clean union');
pause;
