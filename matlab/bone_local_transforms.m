function locals = bone_local_transforms(mesh, globals)
%BONE_LOCAL_TRANSFORMS Compute local transforms for each bone
%
% Inputs
%  mesh: complete mesh object
%  globals: 18 by 1 array of 4 by 4 double returned by bone_axes
%
% Outputs
% locals: 18 by 1 array of 4 by 4 double
%

    locals = cell(size(globals));
    for i = 1 : length(globals)
        p = mesh.bones{i}.parent;
        if p > 0
            locals{i} = inv(globals{p}) * globals{i};
        else
            locals{i} = globals{i};
        end
    end

end

