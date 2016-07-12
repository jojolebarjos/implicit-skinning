function [distance, gradient] = hrbf_apply(centers, coefficients, points)
%HRBF_APPLY Compute HRBF distance value for specified point
%
% Inputs
%  centers: m by 3 double
%  coefficients: m by 4 double
%  points: n by 3 double
%
% Outputs
%  distance: n by 1 double
%  gradient: n by 3 double
%
    
    % Initialize
    m = size(centers, 1);
    n = size(points, 1);
    distance = zeros(n, 1);
    gradient = zeros(n, 3);
    
    % For each center...
    for i = 1 : m
        delta = points - repmat(centers(i, :), n, 1);
        norm = sqrt(sum(delta .* delta, 2));
        
        % Compute distance
        distance = distance + coefficients(i, 1) * norm .* norm .* norm + ...
            (delta * coefficients(i, 2 : 4)') * 3 .* norm;
        
        % Compute gradient
        if nargout > 1
            norm(abs(norm) < 0.000001) = 1;
            gradient = gradient + ...
                repmat(coefficients(i, 1) * norm, 1, 3) .* delta + ...
                repmat(delta * coefficients(i, 2 : 4)' ./ norm, 1, 3) .* delta + ...
                norm * coefficients(i, 2 : 4);
        end
    end
    gradient = 3 * gradient;

end
