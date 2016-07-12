function transform = constraint_transformation(transform, center, slide_axe)
%CONSTRAINT_TRANSFORMATION Constraint transform to fix defined point (some translation can be allowed along an axis)
    
    % Forbid scaling
    scale = (norm(transform(1, 1 : 3)) + norm(transform(2, 1 : 3)) + norm(transform(3, 1 : 3))) / 3;
    transform(1 : 3, 1 : 3) = transform(1 : 3, 1 : 3) / scale;
    
    % Transform fixed point to get desired delta
    delta = center - apply_matrix(transform, center);
    
    % If axis sliding is allowed, adapt delta
    if nargin > 2 && norm(slide_axe) > 0.1
        slide_axe = normalizerow(slide_axe);
        distance = -delta * slide_axe';
        distance = clamp(distance, 0, 0.05);
        delta = delta + slide_axe * distance;
        % TODO forbid rotation around axis?
    end
    
    % Patch transformation
    transform = [eye(3), delta'; 0, 0, 0, 1] * transform;

end

