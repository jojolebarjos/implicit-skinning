function transforms = bone_transforms(mesh, centers, preserve_distances)
%BONE_TRANSFORMS Compute global transforms for each bone
%
% Inputs
%  mesh: complete mesh object
%  centers: 30 by 1 array of spheres or 30 by 3 double
%  preserve_distances: logical or double (default true)
%
% Outputs
%  transforms: 18 by 1 array of 4 by 4 double
%

    % Define missing arguments
    if nargin < 3
        preserve_distances = true;
    end
    if strcmp(preserve_distances, 'optimize')
        % TODO use optimization to compute ideal transforms
        return;
    end
    if islogical(preserve_distances)
        preserve_distances = double(preserve_distances);
    end
    preserve_distances = clamp(preserve_distances, 0, 1);
    
    % Compute original and new axes
    A = bone_axes(mesh.spheres);
    B = bone_axes(centers);
    
    if preserve_distances > 0
        % TODO need to test this distance preservation
    
        % Compute delta required to preserve distances
        deltas = cell(1, length(mesh.bones));
        for i = 1 : length(mesh.bones)
            parent = mesh.bones{i}.parent;
            if parent <= 0
                deltas{i} = [0; 0; 0];
            else
                distance = norm(A{i}(1 : 3, 4) - A{parent}(1 : 3, 4));
                current = B{i}(1 : 3, 4) - B{parent}(1 : 3, 4);
                expected = current * (distance / sqrt(current' * current));
                deltas{i} = expected - current;
            end
        end

        % Propagate delta to children
        done = false(1, length(mesh.bones));
        while ~all(done)
            for i = 1 : length(mesh.bones)
                parent = mesh.bones{i}.parent;
                if parent <= 0
                    done(i) = true;
                elseif done(parent)
                    deltas{i} = deltas{i} + deltas{parent};
                    done(i) = true;
                end
            end
        end

        % Apply delta
        for i = 1 : length(mesh.bones)
            B{i}(1 : 3, 4) = B{i}(1 : 3, 4) + preserve_distances * deltas{i};
        end
    
    end
    
    % Compute transforms
    transforms = cell(1, length(mesh.bones));
    for i = 1 : length(mesh.bones)
        transforms{i} = B{i} * inv(A{i});
    end

end
