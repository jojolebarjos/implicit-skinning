function g = sdf_clean_union(f1, f2)
%SDF_CLEAN_UNION Compute smooth union of given locally supported functions
%
% Inputs
%  f1: n by 1 double
%  f2: n by 1 double
%
% Outputs
%  g: n by 1 double
%
    
    % Dummy vectorization
    if length(f1) ~= 1
        g = arrayfun(@sdf_clean_union, f1, f2);
        return;
    end

    % Clean union introduced by Pasko and al. for global support
    % g = f1 + f2 - sqrt(f1 .* f1 + f2 .* f2);
    
    % Eq. 9 in A Gradient-Implicit Blend paper by Gourmel and al.
    if (f1 <= 0.5 && (f1 <= 2 * f2 * f2 || f2 <= 2 * f1 * f1)) || ...
            (f1 >= 0.5 && (f1 >= 2 * f2 * f2 || f2 >= 2 * f1 * f1))
        g = max(f1, f2);
    else
        % HACK since fzero might not converge otherwise
        % TODO inaccurate in some corner cases
        seed = max(f1, f2);
        g = fzero(@(C) h(f1, f2, C) - 1, seed);
        if g > 1
            g = 1;
        end
    end

end

% Eq. 10 in A Gradient-Implicit Blend paper
function y = h(f1, f2, C)
    if f1 <= 0.5
        a = f1 - 2 * C * C;
        b = f2 - 2 * C * C;
        c = C  - 2 * C * C;
    else
        a = f1 - sqrt(C / 2);
        b = f2 - sqrt(C / 2);
        c = C  - sqrt(C / 2);
    end
    y = sqrt(a * a + b * b) / c;
end
