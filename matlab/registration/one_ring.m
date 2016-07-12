function neighbors = one_ring(faces, vertex_count)
%ONE_RING Get one-ring neighbors

    if nargin < 2
        vertex_count = max(faces(:));
    end
    
    neighbors = cell(1, vertex_count);
    for f = 1 : size(faces, 1)
        face = faces(f, :);
        neighbors{face(1)} = [neighbors{face(1)}, face(2), face(3)];
        neighbors{face(2)} = [neighbors{face(2)}, face(1), face(3)];
        neighbors{face(3)} = [neighbors{face(3)}, face(1), face(2)];
    end
    for i = 1 : vertex_count
        neighbors{i} = setdiff(unique(neighbors{i}), i);
    end

end

