function vectors = matrix_apply(matrix, vectors, constant)
%MATRIX_APPLY Apply 3D transformation to vector
%
% Inputs
%  matrix: 4 by 4 double or 3 by 3 double
%  vectors: n by 3 double
%  constant: double (default 1)
%
% Outputs
%  vectors: n by 3 double
%

    % Validate arguments
    if nargin < 3
        constant = 1;
    end
    if length(constant) == 1
        constant = ones(size(vectors, 1), 1) * constant;
    end
    if size(matrix, 1) == 3
        matrix = [matrix, [0; 0; 0]; 0, 0, 0, 1];
    end
    
    % Apply transformation
    w = matrix * [vectors, constant]';
    vectors = w(1:3, :)';

end

