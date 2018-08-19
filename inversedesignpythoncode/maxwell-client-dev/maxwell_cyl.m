%% maxwell_cyl
% Cylinder of constant epsilon/mu within the simulation grid.

%%% Syntax
%
% * |f = maxwell_cyl(center, radius, height)|
%   produces function handle |f| which defines a cylinder centered at |center|
%   of dimensions |radius| and |height|.

%%% Source code
function [cyl_fun] = maxwell_cyl(center, radius, height)

    validateattributes(center, {'double'}, ...
        {'vector', 'numel', 3}, mfilename, 'center');
    validateattributes(radius, {'double'}, ...
        {'scalar', 'positive',}, mfilename, 'radius');
    validateattributes(height, {'double'}, ...
        {'scalar', 'positive',}, mfilename, 'height');

    box_size = [2*radius, 2*radius, height];
        
    function [is_in, bounding_box] = f(x, y, z)
        bounding_box = {center - box_size/2, ...
                        center + box_size/2};
        is_in = ((x-center(1)).^2 + (y-center(2)).^2 < radius^2) & ...
                (abs(z-center(3)) < height/2);
    end

    cyl_fun = @f;
end
