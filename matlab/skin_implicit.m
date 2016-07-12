function mesh = skin_implicit(mesh, transforms, reparam_threshold, varargin)
%SKIN_IMPLICIT Deform mesh using implicit skinning
%
% Inputs
%  mesh: complete mesh
%  transforms: m by 1 array of 4 by 4 double
%  reparam_threshold: double
%  Optional:
%    'linear' or 'dualquat': Specify whether to use linear blend skinning or dual quaternion skinning. Default: 'dualquat'
%    'spheres' or 'hrbf': Specify whether to use spheres or hrbf as the signed distance functions. Default: 'hrbf'
%
% Outputs
%  mesh: complete mesh
%

    % Parse arguments
    use_linear = any(strcmp('linear', varargin));
    use_spheres = any(strcmp('spheres', varargin));
    varargin = varargin(cell2mat(cellfun(@(s) all(~strcmp(s, {'linear', 'dualquat', 'spheres', 'hrbf'})), varargin, 'UniformOutput', false)));
    for i = 1 : length(varargin)
        warning(['Unrecognized option ', varargin{i}]);
    end
    
    % Compute distances at rest
    implicit.iso_values_at_rest = sdf_field(mesh, use_spheres, reparam_threshold);
    
    % Initial guess using standard skinning
    if use_linear
        mesh = skin_linear(mesh, transforms);
    else
        mesh = skin_dualquat(mesh, transforms);
    end
    
    implicit.has_converged = zeros(size(mesh.vertices, 1), 1);
    implicit.beta = zeros(size(mesh.vertices, 1), 1);
    [implicit.iso_values, implicit.gradients] = sdf_field(mesh, use_spheres, reparam_threshold, transforms);

    % Perform the first iteration of surface tracking
    [mesh, implicit] = vertex_projection(mesh, implicit);
    mesh = tangential_relaxation(mesh, implicit);

    mesh2 = mesh;
    max_iterations = 100;
    iterations = 0;
    
    while ~all(implicit.has_converged) && iterations < max_iterations
        
        % Compute the new iso-values and gradients for the vertices that have not converged
        implicit.prev_gradients = implicit.gradients;
        mesh2.vertices = mesh.vertices(~implicit.has_converged, :);
        mesh2.weights = mesh.weights(~implicit.has_converged, :);
        [implicit.iso_values(~implicit.has_converged), implicit.gradients(~implicit.has_converged, :)] = ...
            sdf_field(mesh2, use_spheres, reparam_threshold, transforms);

        % Perform the surface tracking
        [mesh, implicit] = vertex_projection(mesh, implicit);
        mesh = tangential_relaxation(mesh, implicit);
        iterations = iterations + 1;
    end
    if iterations == max_iterations
        disp(['[implicit] Surface tracking stopped after reaching maximum number of iterations: ', num2str(max_iterations)]);
    end

    % Apply laplacian smoothing at the end to avoid sharp edges on the contact regions
    mesh = laplacian_smoothing(mesh, implicit);
end

function [mesh, implicit] = vertex_projection(mesh, implicit)
%VERTEX_PROJECTION Check the convergence of vertices and deform the mesh by doing one iteration of projection
%
% Inputs
%  mesh: complete mesh
%  implicit: additional data used by implicit
%
% Outputs
%  mesh: complete mesh
%  implicit: additional data used by implicit with updated informations about vertex convergence
%

    sigma = 0.35;
    convergence_threshold = 1e-2;
    
    not_converged = ~implicit.has_converged;
    % Gradient change test: check if the direction of the gradient has changed by more than 55 degrees
    if isfield(implicit, 'prev_gradients')
        angle_between_gradients = acos(dot(implicit.gradients(not_converged, :), implicit.prev_gradients(not_converged, :), 2));
        angle_exceed_threshold = (angle_between_gradients > 55 * pi / 180);
        
        implicit.has_converged(not_converged) = angle_exceed_threshold;
        implicit.beta(not_converged) = angle_exceed_threshold;
    end
    
    not_converged = ~implicit.has_converged;
    % Convergence test: check if the vertex has reached its initial iso-value
    absolute_difference = abs(implicit.iso_values(not_converged) - implicit.iso_values_at_rest(not_converged));
    implicit.has_converged(not_converged) = (absolute_difference <= convergence_threshold);

    not_converged = ~implicit.has_converged;
    % Move the vertices along the gradients
    mesh.vertices(not_converged, :) = mesh.vertices(not_converged, :) + sigma .* ...
        repmat(implicit.iso_values(not_converged) - implicit.iso_values_at_rest(not_converged), [1, 3]) .* implicit.gradients(not_converged, :);
end

function mesh = tangential_relaxation(mesh, implicit)
%TANGENTIAL_RELAXATION Deform the mesh by doing one iteration of tangential relaxation
%
% Inputs
%  mesh: complete mesh
%  implicit: additional data used by implicit
%
% Outputs
%  mesh: complete mesh
%
    
    filter = ~implicit.has_converged & ~mesh.implicit.is_on_side;

    vertices_filtered = mesh.vertices(filter, :);
    gradients_filtered = implicit.gradients(filter, :);
    
    [row, col] = find(mesh.implicit.adj_mat(:, filter));

    % Create matrices than contain replicates of the vertices and gradients (if a vertex v has x neighbours, then v will appear x times)
    vertices = vertices_filtered(col, :);
    gradients = gradients_filtered(col, :);
    
    % Get all the neighbours
    neighbours = mesh.vertices(row, :);
    neighbours_projected = project_onto_plane(vertices, gradients, neighbours);
    
    % Weight each neighbour with its mvc coefficient
    neighbours_weighted = repmat(nonzeros(mesh.implicit.mvc(:, filter)), [1, 3]) .* neighbours_projected;
    
    % Sum along the rows, for each vertex
    % accumarray only accepts vectors... So, we pass each column one by one
    centroids(:, 1) = accumarray(col, neighbours_weighted(:, 1));
    centroids(:, 2) = accumarray(col, neighbours_weighted(:, 2));
    centroids(:, 3) = accumarray(col, neighbours_weighted(:, 3));

    % Transform the vertices
    mu = max(zeros(size(centroids, 1), 1), 1 - (abs(implicit.iso_values(filter) - implicit.iso_values_at_rest(filter)) - 1) .^ 4);
    mesh.vertices(filter, :) = repmat((1 - mu), [1, 3]) .* vertices_filtered + repmat(mu, [1, 3]) .* centroids;
end

function vertices = project_onto_plane(points, normals, vertices)
%PROJECT_ONTO_PLANE Project the vertices onto the planes specified by the points and the normals
%
% Inputs
%  points: n by 3 double
%  normals: n by 3 double
%  vertices: n by 3 double
%
% Outputs
%  vertices: n by 3 double
%

    vertices = vertices - repmat(dot(vertices - points, normals, 2), [1, 3]) .* normals;
end

function mesh = laplacian_smoothing(mesh, implicit)
%LAPLACIAN_SMOOTHING Deform the mesh by applying laplacian smoothing on the surfaces detected as contact regions
%
% Inputs
%  mesh: complete mesh
%  implicit: additional data used by implicit
%
% Outputs
%  mesh: complete mesh
%

    [row, col] = find(mesh.implicit.adj_mat);
    degrees = sum(mesh.implicit.adj_mat, 2);
    % Diffusion of Beta
    for j = 1 : 3
        beta_neighbours = implicit.beta(row, :);
        implicit.beta = (implicit.beta + accumarray(col, beta_neighbours)) ./ (degrees + 1);
    end

    neighbours = mesh.vertices(row, :);
    centroids(:, 1) = accumarray(col, neighbours(:, 1)) ./ degrees;
    centroids(:, 2) = accumarray(col, neighbours(:, 2)) ./ degrees;
    centroids(:, 3) = accumarray(col, neighbours(:, 3)) ./ degrees;
    
    mesh.vertices = repmat(1 - implicit.beta, [1, 3]) .* mesh.vertices + repmat(implicit.beta, [1, 3]) .* centroids;
end
