function [d, g] = sphere_triple(c1, r1, c2, r2, c3, r3, p)
%SPHERE_TRIPLE Compute distance and gradient of triple sphere SDF
%
% Inputs
%  c1: 1 by 3 double
%  r1: double
%  c2: 1 by 3 double
%  r2: double
%  c3: 1 by 3 double
%  r3: double
%  p: n by 3 double
%
% Outputs
%  d: n by 1 double
%  g: n by 3 double
%
    
    % TODO check that trisphere is not degenerated (i.e. this is a bisphere)
    
    % Compute tangent plane
    z12 = c1 + (c2 - c1) * r1 / (r1 - r2);
    z13 = c1 + (c3 - c1) * r1 / (r1 - r3);
    l = normalizerow(z12 - z13);
    z = z12 + ((c1 - z12) * l') * l;
    eta = norm(c1 - z);
    sin_beta = r1 / eta;
    nu = sqrt(eta ^ 2 - r1 ^ 2);
    cos_beta = nu / eta;
    f = (c1 - z) / eta;
    h = normalizerow(cross(l, f));
    k = -sin_beta * h + cos_beta * f;
    tn = normalizerow(z + nu * k - c1);
    t1 = c1 + tn * r1;
    t2 = c2 + tn * r2;
    t3 = c3 + tn * r3;
    
    % Compute central plane normal
    cn = normalizerow(cross(c2 - c1, c3 - c1));
    if cn * tn' < 0
        [d, g] = sphere_triple(c2, r2, c3, r3, c1, r1, p);
        return;
    end
    
    % Compute segment normals
    s12 = t2 - t1;
    s23 = t3 - t2;
    s31 = t1 - t3;
    n12 = cross(s12, tn);
    n23 = cross(s23, tn);
    n31 = cross(s31, tn);
    
    % Vectorize code
    n = size(p, 1);
    d = zeros(n, 1);
    g = zeros(n, 3);
    
    % Mirror points that are on wrong side
    m = (p - repmat(c1, n, 1)) * cn';
    mirror = m < 0;
    p(mirror, :) = p(mirror, :) - 2 * m(mirror, 1) * cn;
    
    % Compute delta for each corner
    d1 = p - repmat(t1, n, 1);
    d2 = p - repmat(t2, n, 1);
    d3 = p - repmat(t3, n, 1);
    
    % First segment
    first = d1 * n12' > 0 & ~(d2 * s23' > 0 & d2 * n23' > 0) & ~(d1 * s31' < 0 & d1 * n31' > 0);
    [d(first), g(first, :)] = sphere_double(c1, r1, c2, r2, p(first, :));
    
    % Second segment
    second = ~first & (d2 * n23' > 0) & ~(d3 * s31' > 0 & d3 * n31' > 0);
    [d(second), g(second, :)] = sphere_double(c2, r2, c3, r3, p(second, :));
    
    % Third segment
    third = ~first & ~second & (d3 * n31' > 0);
    [d(third), g(third, :)] = sphere_double(c3, r3, c1, r1, p(third, :));
    
    % In between
    between = ~first & ~second & ~third;
    b = sum(between);
    d(between) = (p(between, :) - repmat(t1, b, 1)) * tn';
    g(between, :) = repmat(tn, b, 1);
    
    % Mirror again if needed
    g(mirror, :) = g(mirror, :) - 2 * (g(mirror, :) * cn') * cn;

end
