function [vertices, faces] = filter_vertices(vertices, faces, keep)
%FILTER_VERTICES Keep only specified vertices and discard unused elements

    % Drop unused faces
    indices = 1 : size(vertices, 1);
    faces = faces(all(ismember(faces, indices(keep)), 2), :);
    
    % Remap vertices
    map = cumsum(keep);
    for f = 1 : size(faces, 1)
        for i = 1 : 3
            faces(f, i) = map(faces(f, i));
        end
    end
    vertices = vertices(keep, :);

end

