function axes = bone_axes(centers)
%BONE_AXES Compute local axes (local to world transform matrix) for each bone
%
% Inputs
%  centers: 30 by 1 array of spheres or 30 by 3 double
%
% Outpus
%  axes: 18 by 1 array of 4 by 4 double
%
 
    % Extract sphere centers if needed
    if iscell(centers)
        tmp = zeros(30, 3);
        for i = 1 : 30
            tmp(i, :) = centers{i}.center;
        end
        centers = tmp;
    end

    % Compute axes for each bone
    axes = cell(1, 18);
    
    % Palm
    root = centers(26, :) * 0.5 + centers(27, :) * 0.5;
    axes{1} = build(root, centers(22, :) - root, centers(25, :) - centers(24, :));
    axes{2} = build(centers(22, :), centers(12, :) - centers(22, :), centers(20, :) - centers(23, :));
    proto = centers(8, :) * 0.5 + centers(23, :) * 0.5;
    axes{15} = build(proto, centers(4, :) - proto, centers(8, :) - centers(4, :));
    
    % Thumb
    % TODO is thumb normal correct?
    n_thumb = -finger_normal(centers(20, :), centers(19, :), centers(17, :), centers(16, :));
    axes{3} = build(centers(20, :), centers(19, :) - centers(20, :), n_thumb);
    axes{4} = build(centers(19, :), centers(18, :) - centers(19, :), n_thumb);
    axes{5} = build(centers(18, :), centers(17, :) - centers(18, :), n_thumb);
    
    % Index
    n_index = -finger_normal(centers(16, :), centers(15, :), centers(13, :), centers(12, :));
    axes{6} = build(centers(16, :), centers(15, :) - centers(16, :), n_index);
    axes{7} = build(centers(15, :), centers(14, :) - centers(15, :), n_index);
    axes{8} = build(centers(14, :), centers(13, :) - centers(14, :), n_index);
    
    % Middle
    n_middle = -finger_normal(centers(12, :), centers(11, :), centers(9, :), centers(8, :));
    axes{9} = build(centers(12, :), centers(11, :) - centers(12, :), n_middle);
    axes{10} = build(centers(11, :), centers(10, :) - centers(11, :), n_middle);
    axes{11} = build(centers(10, :), centers(9, :) - centers(10, :), n_middle);
    
    % Ring
    n_ring = -finger_normal(centers(8, :), centers(7, :), centers(5, :), centers(4, :));
    axes{12} = build(centers(8, :), centers(7, :) - centers(8, :), n_ring);
    axes{13} = build(centers(7, :), centers(6, :) - centers(7, :), n_ring);
    axes{14} = build(centers(6, :), centers(5, :) - centers(6, :), n_ring);
    
    % Pinky
    n_pinky = finger_normal(centers(4, :), centers(1, :), centers(3, :), centers(8, :));
    axes{16} = build(centers(4, :), centers(3, :) - centers(4, :), n_pinky);
    axes{17} = build(centers(3, :), centers(2, :) - centers(3, :), n_pinky);
    axes{18} = build(centers(2, :), centers(1, :) - centers(2, :), n_pinky);
end

function normal = finger_normal(u, v, w, x)
    uv = normalizerow(v - u);
    uw = normalizerow(w - u);
    ux = normalizerow(x - u);
    n = cross(uw, uv);
    l = norm(n);
    if l < 0.001
        n = ux;
    else
        n = n / l;
    end
    if dot(ux, n) < 0
        n = -n;
    end
    normal = l * n + (1 - l) * ux;
end

function axe = build(o, x, y)
    x = normalizerow(x);
    y = normalizerow(y);
    z = normalizerow(cross(x, y));
    y = cross(z, x);
    axe = [x', y', z', o'; 0, 0, 0, 1];
end
