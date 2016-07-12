function m = rotation_matrix(angle, axis, center)
%ROTATION_MATRIX 3D rotation matrix around axis

    % Fill missing arguments
    if nargin < 3
        center = [0, 0, 0];
    end
    if nargin < 2
        axis = [0, 0, 1];
    end
    
    % Translate toward center
    t1 = eye(4);
    t1(1 : 3, 4) = center';
    
    % Rotate around axis
    r = eye(4);
    r(1 : 3, 1 : 3) = axisangle2matrix(axis, angle);
    
    % Translate back
    t2 = eye(4);
    t2(1 : 3, 4) = -center';
    
    % Combine
    m = t1 * r * t2;

end

