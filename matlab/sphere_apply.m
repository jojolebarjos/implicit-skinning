function [d, g] = sphere_apply(spheres, chunks, p)
%SPHERE_APPLY Compute sphere distance for specified points
%
% Inputs
%  spheres: array of sphere
%  chunks: array of integer matrix
%  p: n by 3 double
%
% Outputs
%  d: n by 1 double
%  g: n by 3 double
%

    % Vectorize code
    n = size(p, 1);
    d = inf(n, 1);
    g = zeros(n, 3);
    
    % For each sphere chunk
    for i = 1 : length(chunks)
        
        % Call appropriate function
        switch length(chunks{i})
            case 1
                s = spheres{chunks{i}(1)};
                [di, gi] = sphere_single(s.center, s.radius, p);
            case 2
                s1 = spheres{chunks{i}(1)};
                s2 = spheres{chunks{i}(2)};
                [di, gi] = sphere_double(s1.center, s1.radius, s2.center, s2.radius, p);
            case 3
                s1 = spheres{chunks{i}(1)};
                s2 = spheres{chunks{i}(2)};
                s3 = spheres{chunks{i}(3)};
                [di, gi] = sphere_triple(s1.center, s1.radius, s2.center, s2.radius, s3.center, s3.radius, p);
            otherwise
                warning('invalid chunk');
                di = inf(n, 1);
                gi = zeros(n, 3);
        end
        
        % Keep minimal distance
        better = di < d;
        d(better) = di(better);
        g(better, :) = gi(better, :);
    end

end
