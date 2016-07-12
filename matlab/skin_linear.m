function mesh = skin_linear(mesh, transforms)
%SKIN_LINEAR Deform mesh using linear skinning
%
% Inputs
%  mesh: complete mesh
%  transforms: m by 1 array of 4 by 4 double
%
% Outputs
%  mesh: complete mesh
%

    % The normals are computed with the transpose inverse of the transforms
    transformsNormals = cell(1, length(transforms));
    for i = 1:length(transforms)
        transformsNormals{i} = inv(transforms{i})';
    end
    
    % Update vertices and normals
    for i = 1 : size(mesh.vertices, 1)
        t = interpolate(transforms, mesh.weights(i, :));
        mesh.vertices(i, :) = matrix_apply(t, mesh.vertices(i, :));
        t = interpolate(transformsNormals, mesh.weights(i, :));
        mesh.normals(i, :) = matrix_apply(t, mesh.normals(i, :), 0);
    end
    mesh.normals = normalizerow(mesh.normals);
    
    % Update spheres
    if isfield(mesh, 'spheres')
        for i = 1 : length(mesh.spheres)
            mesh.spheres{i}.center = matrix_apply(transforms{mesh.spheres{i}.bone}, mesh.spheres{i}.center);
        end
    end

end

function transform = interpolate(transforms, weights)
    transform = zeros(4);
    for i = 1 : length(transforms)
        transform = transform + transforms{i} * weights(i);
    end
end
