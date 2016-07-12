function g = sdf_blend(f1, f2, theta)
%SDF_BLEND Compute blend as described in Implicit Skinning
%
% Inputs
%  f1: n by 1 double
%  f2: n by 1 double
%  f3: n by 1 double or double
%
% Outputs
%  g: n by 1 double
%

    % Vectorize code
    g = zeros(size(f1));
    theta = abs(theta);
    t = tan(theta);
    if length(f1) ~= length(t)
        theta = repmat(theta, length(f1), 1);
        t = repmat(t, length(f1), 1);
    end
    
    % Pure union outside trapezoidal region
    union = ...
        theta >= pi / 4 | ...
        f2 <= f1 .* t | ...
        f1 <= f2 .* t | ...
        f2 <= 0.5 * (7 - 5 * t) .* f1 + 7 / 4 * (t - 1) | ...
        f1 <= 0.5 * (7 - 5 * t) .* f2 + 7 / 4 * (t - 1);
    g(union) = max(f1(union), f2(union));
    
    % Solve quadratic equation inside special region
    g(~union) = curve(f1(~union), f2(~union), t(~union));
    
end

function g = curve(f1, f2, t)
    
    % Compute coefficients according to area
    upper = (f1 - 0.5 * t) .^ 2 + (f2 - 0.5 * t) .^ 2 > (0.5 - 0.5 * t) .^ 2 & (f1 > t * 0.5 | f2 > t * 0.5);
    a = zeros(size(f1));
    b = -t;
    a(upper) = -7 / 4 * (t(upper) - 1);
    b(upper) = -0.5 * (7 - 5 * t(upper));
    
    % Solve quadratic equation
    u = b .* b - 2 * b - 1;
    v = 2 * b .* (f1 + f2 + a) - 2 * a;
    w = (f1 + a) .* (f1 + a) + (f2 + a) .* (f2 + a) - a .* a;
    singular = abs(u) < 0.0001;
    g = zeros(size(f1));
    g(singular) = -w(singular) ./ v(singular);
    g(~singular) = (-v(~singular) - sqrt(v(~singular) .* v(~singular) - 4 * u(~singular) .* w(~singular))) ./ (2 * u(~singular));
    
end
