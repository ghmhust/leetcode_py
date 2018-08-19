%% maxwell_cyl_smooth
% Cylinder of constant epsilon/mu within the simulation grid with smoothed boundaries.

%%% Syntax
%
% * |f = maxwell_cyl(center, radius, height)|
%   produces function handle |f| which defines a cylinder centered at |center|
%   of dimensions |radius| and |height|.
%
% * |f = maxwell_cyl_smooth(..., 'smooth_dist', d)|
%   smoothes the edges of the cylinder over a distance |d| from its boundary.
%

%%% Source code
function [cyl_fun] = maxwell_cyl_smooth(center, radius, height, varargin)


        %
        % Validate and parse inputs.
        %

    validateattributes(center, {'double'}, {'numel', 3, 'nonnan', 'finite'}, ...
                        mfilename, 'center');

    validateattributes(radius, {'double'}, {'scalar', 'positive'}, ...
                        mfilename, 'radius');

    validateattributes(height, {'double'}, {'scalar', 'positive'}, ...
                        mfilename, 'height');

    % Optional parameters.
    options = my_parse_options(struct('smooth_dist', 1), ...
                                varargin, mfilename);
    validateattributes(options.smooth_dist, {'numeric'}, ...
        {'positive', 'scalar'}, mfilename, 'smooth_dist');


        %
        % Create function handle.
        %

    box_size = [2*radius, 2*radius, height] + 2*options.smooth_dist;
        
    function [is_in, bounding_box] = f(x, y, z)
        bounding_box = {center - box_size/2, ...
                        center + box_size/2};
        dist_r = (radius - sqrt((x-center(1)).^2 + (y-center(2)).^2));
        dist_z = (height/2 - abs(z-center(3)));
        is_in = my_val_clamp(dist_r, options.smooth_dist) .* ...
                my_val_clamp(dist_z, options.smooth_dist);
    end

    cyl_fun = @f;
end


