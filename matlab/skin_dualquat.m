function mesh = skin_dualquat(mesh, transforms)
%SKIN_DUALQUAT Deform mesh using dual quaternion skinning
%
% Inputs
%  mesh: complete mesh
%  transforms: m by 1 array of 4 by 4 double
%
% Outputs
%  mesh: complete mesh
%

    % Initialize
    bones_count = length(transforms);
    Q = zeros(bones_count, 4);
    T = zeros(bones_count, 3);
    
    % For the formulas on how to convert a rotation matrix to a quaternion, see:
    % http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
    % Additionally, it is more robust to compare the traces:
    % http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/ethan.htm
    % The rotation part of a transformation t should have det(t) = 1 and be orthogonal: t' = inv(t).
    for i = 1 : bones_count
        rotation = transforms{i}(1 : 3, 1 : 3);
        translation = transforms{i}(1 : 3, 4);
        
        trace1 = trace(rotation);
        trace2 = rotation(1, 1) - rotation(2, 2) - rotation(3, 3);
        trace3 = rotation(2, 2) - rotation(1, 1) - rotation(3, 3);
        trace4 = rotation(3, 3) - rotation(1, 1) - rotation(2, 2);
        
        [~, index] = max([trace1, trace2, trace3, trace4]);
        switch index
            case 1
                S = 2 * sqrt(1 + trace1);
                qr = 0.25 * S;
                qi = (rotation(3, 2) - rotation(2, 3)) / S;
                qj = (rotation(1, 3) - rotation(3, 1)) / S;
                qk = (rotation(2, 1) - rotation(1, 2)) / S;
            case 2
                S = 2 * sqrt(1 + trace2);
                qr = (rotation(3, 2) - rotation(2, 3)) / S;
                qi = 0.25 * S;
                qj = (rotation(1, 2) + rotation(2, 1)) / S;
                qk = (rotation(1, 3) + rotation(3, 1)) / S;
            case 3
                S = 2 * sqrt(1 + trace3);
                qr = (rotation(1, 3) - rotation(3, 1)) / S;
                qi = (rotation(1, 2) + rotation(2, 1)) / S;
                qj = 0.25 * S;
                qk = (rotation(2, 3) + rotation(3, 2)) / S;
            case 4
                S = 2 * sqrt(1 + trace4);
                qr = (rotation(2, 1) - rotation(1, 2)) / S;
                qi = (rotation(1, 3) + rotation(3, 1)) / S;
                qj = (rotation(2, 3) + rotation(3, 2)) / S;
                qk = 0.25 * S;
        end
        quat = [qr, qi, qj, qk];
        
        parent = mesh.bones{i}.parent;
        % If the dot product between this quaternion and the one of the parent is negative, change the sign of this quaternion
        if parent ~= 0
            % TODO: Is it safe to assume that parent < i? If not, we need to re-order/do multiple pass
            if parent > i
                error('Unhandled topology: The root should be computed before the leaves');
            end
            if dot(Q(parent, :), quat) < 0
                quat = -quat;
            end
        end
        
        Q(i, :) = quat;
        T(i, :) = translation';
    end
    
    DQ = quattrans2udq(Q, T);
    
    % The normals can be computed by evaluating only the non-dual part of the dual quaternion, see here for more:
    % http://dev.theomader.com/dual-quaternion-skinning/
    DQ_normals = DQ;
    for i = 1:bones_count
        % Get rid of the dual part
        DQ_normals(2, :, i) = zeros(1, 1, 4);
    end
    
    % Update spheres: As the spheres are only attached to one bone, dual quat is equivalent to linear.
    % So, we can just apply directly the linear transform matrices.
    if isfield(mesh, 'spheres')
        for i = 1 : length(mesh.spheres)
            mesh.spheres{i}.center = matrix_apply(transforms{mesh.spheres{i}.bone}, mesh.spheres{i}.center);
        end
    end
    
    mesh.vertices = dualquatlbs(mesh.vertices, DQ, mesh.weights);
    mesh.normals = dualquatlbs(mesh.normals, DQ_normals, mesh.weights);
    
end
