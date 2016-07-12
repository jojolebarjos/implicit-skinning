function g = sdf_union(f1, f2)
%SDF_UNION Compute union of given locally supported functions
%
% Inputs
%  f1: n by 1 double
%  f2: n by 1 double
%
% Outputs
%  g: n by 1 double
%

    % Ricci's operator is the maximum
    g = max(f1, f2);

end
