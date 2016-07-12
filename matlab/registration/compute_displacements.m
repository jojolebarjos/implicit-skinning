function displacements = compute_displacements(displacements, pairs, mesh_vertices, mesh_neighbors, points_vertices, points_normals, weight_smoothness, weight_closeness)
%COMPUTE_DISPLACEMENTS Compute incremental vertex displacements

    % Define missing parameters
    if nargin < 7
        weight_smoothness = 1;
    end
    if nargin < 8
        weight_closeness = 1e-2;
    end
    N = size(pairs, 1);
    M = size(mesh_vertices, 1);
    
    % Build distance matrix and right-hand term
    I = zeros(1, 3 * N);
    J = zeros(1, 3 * N);
    V = zeros(1, 3 * N);
    b = zeros(N, 1);
    for i = 1 : N
        j = pairs(i, 1);
        p = mesh_vertices(j, :);
        q = points_vertices(pairs(i, 2), :);
        n = points_normals(pairs(i, 2), :);
        
        I(3 * i - 2) = i;
        J(3 * i - 2) = 3 * j - 2;
        V(3 * i - 2) = -n(1);
        
        I(3 * i - 1) = i;
        J(3 * i - 1) = 3 * j - 1;
        V(3 * i - 1) = -n(2);
        
        I(3 * i) = i;
        J(3 * i) = 3 * j;
        V(3 * i) = -n(3);
        
        b(i) = n * (q - p)';
        
    end
    J_distance = sparse(I, J, V, N, 3 * M);
    
    % Build closeness matrix (identity)
    J_closeness = speye(3 * M, 3 * M);
    
    % Build laplacian smoothing matrix
    I = [];
    J = [];
    V = [];
    for i = 1 : N
        
        % Add coefficient for each neighbor
        for j = mesh_neighbors{i}
            I = [I, 3 * i - 2, 3 * i - 1, 3 * i];
            J = [J, 3 * j - 2, 3 * j - 1, 3 * j];
            V = [V, 1, 1, 1];
        end
        
        % Add coefficient for vertex
        I = [I, 3 * i - 2, 3 * i - 1, 3 * i];
        J = [J, 3 * i - 2, 3 * i - 1, 3 * i];
        valence = length(mesh_neighbors{i});
        V = [V, -valence, -valence, -valence];
        
    end
    J_smoothness = sparse(I, J, V, 3 * M, 3 * M);
    
    % Create and solve linear system
    JtJ = (J_distance' * J_distance) + weight_closeness * J_closeness + weight_smoothness * (J_smoothness' * J_smoothness);
    Jtb = J_distance' * (J_distance * reshape(displacements', 3 * M, 1) - b);
    displacements = reshape(JtJ \ Jtb, M, 3);

end
