function holes = hole_find(faces)
%HOLE_FIND Find all holes in given mesh (i.e. sequence of boundary vertices)

    % Find boundary edges (assuming CCW faces)
    edges = [ ...
        faces(:, [1, 2]); ...
        faces(:, [2, 3]); ...
        faces(:, [3, 1]) ...
    ];
    borders = [];
    for i = 1 : size(edges, 1)
        if ~any(edges(:, 1) == edges(i, 2) & edges(:, 2) == edges(i, 1))
            borders = [borders; edges(i, [2, 1])];
        end
    end

    % Extract hole sequences
    holes = {};
    while ~isempty(borders)
        hole = borders(1, :);
        borders = borders(2 : end, :);
        while hole(1) ~= hole(end)
            c = find(borders(:, 1) == hole(end));
            hole = [hole, borders(c, 2)];
            borders(c, :) = [];
        end
        hole = hole(2 : end)';
        holes = [holes, hole];
    end

end
