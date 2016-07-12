function theta = sdf_blend_theta_curve(alpha, alpha0, alpha1, alpha2, theta0, theta1, theta2, omega0, omega1)
%SDF_BLEND_THETA_CURVE Map alpha to theta as described in Implicit Skinning
%
% Inputs
%  alpha: n by 1 double
%  alpha0: double
%  alpha1: double
%  alpha2: double
%  theta0: double
%  theta1: double
%  theta2: double
%  omega0: double
%  omega1: double (default omega0)
%
% Outputs
%  theta: n by 1 double
%

    % Fill missing arguments
    if nargin == 2
        alpha1 = alpha0.alpha1;
        alpha2 = alpha0.alpha2;
        theta0 = alpha0.theta0;
        theta1 = alpha0.theta1;
        theta2 = alpha0.theta2;
        omega0 = alpha0.omega0;
        omega1 = alpha0.omega1;
        alpha0 = alpha0.alpha0;
    elseif nargin < 9
        omega1 = omega0;
    end
    
    % Split input
    a = alpha <= alpha0;
    b = alpha > alpha0 & alpha <= alpha1;
    c = alpha > alpha1 & alpha < alpha2;
    d = alpha >= alpha2;
    
    % Piecewise interpolation
    theta = zeros(size(alpha));
    theta(a) = theta0;
    theta(b) = kappa((alpha(b) - alpha1) / (alpha0 - alpha1)) .^ omega0 * (theta0 - theta1) + theta1;
    theta(c) = kappa((alpha(c) - alpha1) / (alpha2 - alpha1)) .^ omega1 * (theta2 - theta1) + theta1;
    theta(d) = theta2;
    
end

function y = kappa(x)
    x = max(min(x, 1), 0);
    x(x == 0) = 0; % Negative zeros are evil!!!
    y = 1 - exp(1 - 1 ./ (1 - exp(1 - 1 ./ x)));
end
