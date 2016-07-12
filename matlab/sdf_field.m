function [field, gradients] = sdf_field(mesh, use_spheres, threshold, transforms)
%SDF_FIELD Compute the isovalues and the gradients at the vertices positions
%
% Inputs
%  mesh: complete mesh
%  use_sphere: logical, Specify whether to use spheres or hrbf
%  reparam_threshold: double
%  Optional:
%    transforms: m by 1 array of 4 by 4 double, Unused when using spheres, Set to identity if not provided
%
% Outputs
%  field: n by 1 double
%  gradients: n by 3 double
%

    % Select SDF function
    if use_spheres
        sdf_func = @(bone, vertices) sphere_apply(mesh.spheres, mesh.bones{bone}.blocks, vertices);
    else
        % If the transforms are not provided, set them to the identity
        if nargin < 4
            transforms = cellfun(@double, repmat({eye(4)}, 1, length(mesh.bones)), 'Un', 0);
        end
        sdf_func = @hrbf_field;
    end
    
    function [dist, grad] = hrbf_field(bone, vertices)
        % Compute the hrbf distance and gradients in the reference frame of the rest pose
        vertices_transformed = matrix_apply(inv(transforms{bone}), vertices);
        [dist, grad] = hrbf_apply(mesh.bones{bone}.hrbf.centers, mesh.bones{bone}.hrbf.coefficients, vertices_transformed);
        % Change the reference frame of the gradients back
        grad = matrix_apply(transforms{bone}, grad, 0);
    end

    if ~isfield(mesh.implicit, 'alpha_to_theta')
        % If alpha to theta parameters are not provided, set them to default values
        mesh.implicit.alpha_to_theta.alpha0 = pi / 8;
        mesh.implicit.alpha_to_theta.alpha1 = pi / 2;
        mesh.implicit.alpha_to_theta.alpha2 = pi / 2;
        mesh.implicit.alpha_to_theta.theta0 = pi / 4;
        mesh.implicit.alpha_to_theta.theta1 = pi / 10;
        mesh.implicit.alpha_to_theta.theta2 = pi / 10;
        mesh.implicit.alpha_to_theta.omega0 = 0.5;
        mesh.implicit.alpha_to_theta.omega1 = 1;
    end
    
    % For each bone, compute the signed distance function at the points
    sdf = cell(length(mesh.bones), 1);
    for child = 1 : length(mesh.bones)
        parent = mesh.bones{child}.parent;
        distances = zeros(size(mesh.vertices, 1), 4);
        
        % Project the vertices of the child and its parent onto this bone
        non_zero_child = mesh.weights(:, child) ~= 0;
        if parent ~= 0
            non_zero_parent = mesh.weights(:, parent) ~= 0;
        else
            non_zero_parent = false(size(mesh.vertices, 1), 1);
        end
        non_zero_both = non_zero_child | non_zero_parent;
        
        [dist, grad] = sdf_func(child, mesh.vertices(non_zero_both, :));

        % Reparametrize the distances to have a local support in [0, 1] and normalize the gradients
        distances(non_zero_both, :) = [sdf_reparam(dist, threshold), normalizerow(grad)];
        sdf{child} = distances;
        
        if parent ~= 0
            % Project the vertices of the child onto the parent
            [dist, grad] = sdf_func(parent, mesh.vertices(non_zero_child, :));
            sdf{parent}(non_zero_child, :) = [sdf_reparam(dist, threshold), normalizerow(grad)];
        end
    end
    
    % For each vertex, we compute the blending between pairs of bones, then take the union (max)
    [field, gradients] = blend(mesh, sdf, sdf_func, threshold);
    
    if nargout > 1 && any(any(isnan(gradients)))
        % If the reparametrization threshold is too small, then some vertices have a distance of 0 and
        % the gradients are undefined.
        error('Found gradients with value of NaN. Maybe you should retry with a higher threshold value.');
    end
    
end

function [fields, gradients] = blend(mesh, sdf, sdf_func, threshold)

    % Vectorize code
    fields = zeros(size(mesh.vertices, 1), 1);
    gradients = zeros(size(mesh.vertices, 1), 3);
    f = fields;
    g = gradients;
    
    % For each bone (except root)
    for child = 1 : length(mesh.bones)
        parent = mesh.bones{child}.parent;
        if parent == 0
            continue;
        end

        % For each non-zero weighted vertex
        non_zero_child = mesh.weights(:, child) ~= 0;
        non_zero_parent = mesh.weights(:, parent) ~= 0;
        non_zero_both = non_zero_child | non_zero_parent;
        
        if strcmp(mesh.bones{child}.blend, 'union')
            [f(non_zero_both), g(non_zero_both, :)] = union_op( ...
                sdf{child}(non_zero_both, 1), ...
                sdf{child}(non_zero_both, 2 : 4), ...
                sdf{parent}(non_zero_both, 1), ...
                sdf{parent}(non_zero_both, 2 : 4) ...
            );
        else
%             [f_union, g_union] = union_op( ...
%                 sdf{child}(non_zero_both, 1), ...
%                 sdf{child}(non_zero_both, 2 : 4), ...
%                 sdf{parent}(non_zero_both, 1), ...
%                 sdf{parent}(non_zero_both, 2 : 4) ...
%             );
            [f_blend, g_blend] = blend_op( ...
                sdf{child}(non_zero_both, 1), ...
                sdf{child}(non_zero_both, 2 : 4), ...
                sdf{parent}(non_zero_both, 1), ...
                sdf{parent}(non_zero_both, 2 : 4), ...
            child, parent, mesh.vertices(non_zero_both, :), sdf_func, threshold, mesh.implicit.alpha_to_theta);

            % Using the blend weights produces nice, smooth fields, but the other parameters have almost no effects
%             f(non_zero_both, :) = ... 
%                 (1 - mesh.implicit.blend_weight(non_zero_both, :)) .* f_union ...
%                 + mesh.implicit.blend_weight(non_zero_both, :) .* f_blend;
%             g(non_zero_both, :) = slerp(g_union, g_blend, mesh.implicit.blend_weight(non_zero_both, :));
            f(non_zero_both, :) = f_blend;
            g(non_zero_both, :) = g_blend;
            
        end
        
        % Merge with current result
        [fields, gradients] = union_op(fields, gradients, f, g);
        
    end
    
end

function [f, g] = union_op(f1, g1, f2, g2)
    f = f1;
    g = g1;
    k = f2 > f1;
    f(k) = f2(k);
    g(k, :) = g2(k, :);
end

function [f, g] = blend_op(f1, g1, f2, g2, b1, b2, v, sdf_func, threshold, alpha_to_theta)
    
    % Compute distances
    alpha = acos(dot(g1, g2, 2));
    
    % Tweaked Vaillant Contact (It has a shape like figure 7a)
    theta = sdf_blend_theta_curve(alpha, alpha_to_theta);
    f = sdf_blend(f1, f2, theta);
    
    % Get close neighbors to compute numerical derivatives
    neighbors = kron(v, ones(6, 1)) + repmat([ ...
        -1,  0,  0; ...
         1,  0,  0; ...
         0, -1,  0; ...
         0,  1,  0; ...
         0,  0, -1; ...
         0,  0,  1  ...
    ] * 1e-5, length(f1), 1);
    
    % Compute neighbors distances
    % NOTE: we assume same theta for all neighbors
    d1 = sdf_reparam(sdf_func(b1, neighbors), threshold);
    d2 = sdf_reparam(sdf_func(b2, neighbors), threshold);
    blended = sdf_blend(d1, d2, kron(theta, ones(6, 1)));
    
    % Compute numerical gradients
    g = reshape(blended, 6, length(f1))' * [1, 0, 0; -1, 0, 0; 0, 1, 0; 0, -1, 0; 0, 0, 1; 0, 0, -1];
    g = normalizerow(g);
    
end
