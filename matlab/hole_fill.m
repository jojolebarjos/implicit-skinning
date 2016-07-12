function [vertices, faces] = hole_fill(vertices, faces, hole, inflate)
%HOLE_FILL Fill specified hole (CCW face vertex order)

    % Fill missing arguments
    if nargin < 4
        inflate = true;
    end
    
    % Fit initial patch
    [PV, PF] = create_patch(vertices, hole);
    
    if inflate
        
        % Estimate normals
        PN = diffuse_normals(vertices, faces, PV, PF, hole);

        % Adjust vertex locations
        % TODO this does not work well :'(
        PV = adjust_vertices(PV, PF, 1 : length(hole), PN);
    
    end
    
    % Merge geometry
    map = [hole; (1 : size(PV, 1) - length(hole))' + size(vertices, 1)];
    faces = [faces; map(PF)];
    vertices = [vertices; PV(length(hole) + 1 : end, :)];
    
end

function [PV, PF] = create_patch(V, H)
    h = length(H);
    
    % Generate angles and radii
    T = [h; (floor(h / 6) * 6 - 6 : -6 : 6)'; 1];
    R = [];
    A = [];
    for i = 1 : length(T)
        R = [R; ones(T(i), 1) * (length(T) - i)];
        A = [A; (1 : T(i))' / T(i) * 2 * pi];
    end
    
    % Create circular patch with average vertex valence of 6
    PV = repmat(R, 1, 2) .* [cos(A), sin(A)];
    PF = delaunay(PV);
    PV = [PV, zeros(size(PV, 1), 1)];
    %show();
    
    % Get patch segments
    PS = [ ...
        PF(:, [1, 2]); ...
        PF(:, [1, 3]); ...
        PF(:, [2, 3]) ...
    ];
    PS = sort(PS, 2);
    PS = unique(PS, 'rows');
    
    % Compute laplacian matrix
    L = [ ...
        PS(:, 1), PS(:, 1), ones(size(PS, 1), 1); ...
        PS(:, 2), PS(:, 2), ones(size(PS, 1), 1); ...
        PS(:, 1), PS(:, 2), -ones(size(PS, 1), 1); ...
        PS(:, 2), PS(:, 1), -ones(size(PS, 1), 1) ...
    ];
    L = sparse(L(:, 1), L(:, 2), L(:, 3));
    
    % Use quadratic program to find initial guess
    Q = kron(L, speye(3));
    E = speye(3 * h, 3 * size(PV, 1));
    d = reshape(V(H, :)', 3 * h, 1);
    x = [Q, E'; E, sparse(3 * h, 3 * h)] \ [zeros(3 * size(PV, 1), 1); d];
    PV = reshape(x(1 : 3 * size(PV, 1)), 3, size(PV, 1))';
    %show();
    
    % Debug function
%     function show()
%         trimesh(PF, PV(:, 1), PV(:, 2), PV(:, 3));
%         hold on;
%         N = normalizerow(normals(PV, PF));
%         C = (PV(PF(:, 1), :) + PV(PF(:, 2), :) + PV(PF(:, 3), :)) / 3;
%         quiver3(C(:, 1), C(:, 2), C(:, 3), N(:, 1), N(:, 2), N(:, 3))
%         hold off;
%         view([-90, 0]);
%         camlight;
%         view([90, 0]);
%         camlight;
%         axis equal;
%         grid off;
%         lighting gouraud;
%         axis off;
%         pause;
%     end

end

function PN = diffuse_normals(V, F, PV, PF, H)

    % Get patch segments
    PS = [ ...
        PF(:, [1, 2]); ...
        PF(:, [1, 3]); ...
        PF(:, [2, 3]) ...
    ];
    PS = sort(PS, 2);
    PS = unique(PS, 'rows');
    
    % Compute laplacian matrix
    L = [ ...
        PS(:, 1), PS(:, 1), ones(size(PS, 1), 1); ...
        PS(:, 2), PS(:, 2), ones(size(PS, 1), 1); ...
        PS(:, 1), PS(:, 2), -ones(size(PS, 1), 1); ...
        PS(:, 2), PS(:, 1), -ones(size(PS, 1), 1) ...
    ];
    L = sparse(L(:, 1), L(:, 2), L(:, 3));
    
    % Compute smoothing matrix
    S = abs(L);
    S = S' + S;
    
    % Get mesh normals
    N = normalizerow(per_vertex_normals(V, F));
    
    % Get initial guess for patch normals
    PN = [N(H, :); zeros(size(PV, 1) - length(H), 3)];
    TN = normalizerow(per_vertex_normals(PV, PF));
    PN(end, :) = TN(end, :);
    
    % Estimate desired patch vertex normals by diffusion
    IN = PN;
    II = [1 : length(H), size(PV, 1)];
    for i = 1 : length(H)
        PN = S * PN;
        l = sqrt(sum(PN .^ 2, 2));
        l(l < 0.0001) = 1;
        PN = PN ./ repmat(l, 1, 3);
        PN(II, :) = IN(II, :);
        %show();
    end
    
    % Compute desired per face normals
    PN = normalizerow( ...
        PN(PF(:, 1), :) + ...
        PN(PF(:, 2), :) + ...
        PN(PF(:, 3), :) ...
    );

    % Debug function
%     function show()
%         trimesh(PF, PV(:, 1), PV(:, 2), PV(:, 3));
%         hold on;
%         quiver3(PV(:, 1), PV(:, 2), PV(:, 3), PN(:, 1), PN(:, 2), PN(:, 3))
%         hold off;
%         view([-90, 0]);
%         camlight;
%         view([90, 0]);
%         camlight;
%         axis equal;
%         grid off;
%         lighting gouraud;
%         axis off;
%         pause;
%     end

end

function V = adjust_vertices(V, F, H, N)
    f = size(F, 1);
    v = size(V, 1);
    h = length(H);
    
    % Estimate edge size
    S = [ ...
        F(:, [1, 2]); ...
        F(:, [1, 3]); ...
        F(:, [2, 3]) ...
    ];
    r = pi / 2 * mean(sqrt(sum((V(S(:, 1), :) - V(S(:, 2), :)) .^ 2, 2)));
    
    % Compute orthonormal basis {e1, e2, n} for each face
    E2 = rand(f, 3) - 0.5;
    E1 = normalizerow(cross(E2, N, 2));
    E2 = normalizerow(cross(N, E1, 2));
    
    % Get subvertices (face vertices)
    S = reshape(F', 3 * f, 1); 
    P = V(S, :);
    
    % Repeat multiple times
    for n = 1 : 3
    
        % Project subvertices on oriented circle
        C = (V(F(:, 1), :) + V(F(:, 2), :) + V(F(:, 3), :)) / 3;
        D = P - kron(C, [1; 1; 1]);
        PC = [dot(D, kron(E1, [1; 1; 1]), 2), dot(D, kron(E2, [1; 1; 1]), 2)];
        R = ( ...
            sqrt(sum((D(:, 2) - D(:, 1)) .^ 2, 2)) + ...
            sqrt(sum((D(:, 3) - D(:, 2)) .^ 2, 2)) + ...
            sqrt(sum((D(:, 1) - D(:, 3)) .^ 2, 2)) + ...
            r ...
        ) / 4;
        PC = normalizerow(PC) .* repmat(R, 1, 2);
        D = repmat(PC(:, 1), 1, 3) .* kron(E1, [1; 1; 1]) + repmat(PC(:, 2), 1, 3) .* kron(E2, [1; 1; 1]);
        P = D + kron(C, [1; 1; 1]);
        show();

        % Use quadratic program to estimate best centers and vertices
        A = [];
        for k = 1 : 3 * f
            i = ceil(k / 3);
            j = S(k);
            I = (1 : 3) + 3 * k - 3;
            J = [3 * i - 2, 3 * i - 1, 3 * i, 3 * f + 3 * j - 2, 3 * f + 3 * j - 1, 3 * f + 3 * j, 3 * f + 3 * v + 1];
            [J, I] = meshgrid(J, I);
            K = [eye(3), -eye(3), D(k, :)'];
            A = [A; I(:), J(:), K(:)];
        end
        A = sparse(A(:, 1), A(:, 2), A(:, 3), 9 * f, 3 * f + 3 * v + 1);
        E = [zeros(3 * h, 3 * f), kron(sparse(1 : h, H, ones(1, h), h, v), speye(3)), zeros(3 * h, 1); zeros(1, 3 * f + 3 * v), 1];
        Q = [A' * A, E'; E, zeros(3 * h + 1, 3 * h + 1)];
        d = [zeros(3 * f + 3 * v + 1, 1); reshape(V(H, :)', 3 * h, 1); 1];
        x = Q \ d;

        % Update centers
        C = reshape(x(1 : 3 * f), 3, f)';
        P = D + kron(C, [1; 1; 1]);
        show();

        % Update vertices
        V = reshape(x(3 * f + 1 : 3 * f + 3 * v), 3, v)';
        P = V(S, :);
        show();
    
    end
    
    % Debug function
    function show()
        trimesh(reshape(1 : 3 * f, 3, f)', P(:, 1), P(:, 2), P(:, 3));
        hold on;
        quiver3(C(:, 1), C(:, 2), C(:, 3), N(:, 1), N(:, 2), N(:, 3));
        hold off;
        view([-90, 0]);
        camlight;
        view([90, 0]);
        camlight;
        axis equal;
        grid off;
        lighting gouraud;
        axis off;
        pause;
    end
    
end
