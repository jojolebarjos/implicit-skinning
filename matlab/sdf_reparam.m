function distance = sdf_reparam(distance, threshold)
%SDF_REPARAM Reparametrize 0-based distance value to 0.5-based SDF, as explained in the paper Implicit Skinning
%
% Inputs
%  distance: n by 1 double
%  threshold: double
%
% Outputs
%  distance: n by 1 double
%

    lowerThan = distance < -threshold;
    higherThan = distance > threshold;
    inRange = not(lowerThan | higherThan);
    
    distance(lowerThan) = 1;
    distance(higherThan) = 0;
    
    distance(inRange) = distance(inRange) / threshold;
    distance(inRange) = -3 / 16 * distance(inRange) .^ 5 ...
        + 5 / 8 * distance(inRange) .^ 3 ...
        - 15 / 16 * distance(inRange) + 0.5;

end
