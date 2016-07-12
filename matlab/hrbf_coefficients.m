function coefficients = hrbf_coefficients(centers, points, normals)
%HRBF_COEFF Compute HRBF coefficients for provided centers
%
% Inputs
%  centers: m by 3 double
%  points: n by 3 double
%  normals: n by 3 double
%
% Outputs
%  coefficients: m by 4 double
%

    N = size(points, 1);
    M = size(centers, 1);

    % Build linear system
    A = zeros(4 * N, 4 * M);
    B = zeros(4 * N, 1);
    a = zeros(4, 4);
    b = zeros(4, 1);
    for n = 1 : N
        
        % Fill right-hand side
        b(1) = 0;
        b(2 : 4) = normals(n, :);
        B(4 * n - 3 : 4 * n) = b;
        
        % Fill left-hand side
        for m = 1 : M
            
            x = points(n, :);
            v = centers(m, :);
            delta = x - v;
            norm = sqrt(delta * delta');
            
            if norm > 0.00001
            
                phi = norm * norm * norm;
                dphi = 3 * norm * norm;
                ddphi = 6 * norm;

                a(1, 1) = phi;
                a(1, 2 : 4) = delta * dphi / norm;
                a(2 : 4, 1) = delta * dphi / norm;
                a(2 : 4, 2 : 4) = (ddphi - dphi / norm) / (norm * norm) * (delta' * delta) + eye(3) * dphi / norm;
                
                A(4 * n - 3 : 4 * n, 4 * m - 3 : 4 * m) = a;
                
            end
        end
    end
    
    % Solve linear system
    coefficients = reshape(A \ B, [4, M])';

end
