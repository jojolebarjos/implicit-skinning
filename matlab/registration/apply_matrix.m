function v = apply_matrix(m, v, constant)
%APPLY_MATRIX apply 3D transformation to vector

    if nargin < 3
        constant = 1;
    end
    if length(constant) == 1
        constant = ones(size(v, 1), 1) * constant;
    end
    if size(m, 1) == 3
        m = [m, [0; 0; 0]; 0, 0, 0, 1];
    end
    
    w = m * [v, constant]';
    v = w(1:3, :)';

end

