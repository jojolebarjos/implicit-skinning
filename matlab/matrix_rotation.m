function matrix = matrix_rotation(angle, axis, center)
%MATRIX_ROTATION 3D rotation matrix around axis
%
% Inputs
%  angle: double
%  axis: 1 by 3 double (default [0, 0, 1])
%  center: 1 by 3 double (default [0, 0, 0])
%
% Outputs
%  matrix: 4 by 4 double
%

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
    matrix = t1 * r * t2;

end

