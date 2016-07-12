function [d, g] = sphere_single(c, r, p)
%SPHERE_SINGLE Compute distance and gradient of single sphere SDF
%
% Inputs
%  c: 1 by 3 double
%  r: double
%  p: n by 3 double
%
% Outputs
%  d: n by 1 double
%  g: n by 3 double
%

    n = size(p, 1);
    delta = p - repmat(c, n, 1);
    norm = sqrt(sum(delta .* delta, 2));
    d = norm - r;
    g = delta ./ repmat(norm, 1, 3);

end

