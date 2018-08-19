%% maxwell_box_smooth
% Box of constant epsilon/mu within the simulation grid with smoothed boundaries.

%%% Syntax
%
% * |fun = maxwell_box_smooth(center, box_size)| 
%   returns function handle |fun| which describes a rectangular prism 
%   centered at |center| and of size |box_size|.
%
% * |fun = maxwell_box_smooth(..., 'smooth_dist', d)| 
%   smoothes the edges of the box over a distance |d|.
%


%%% Source code
function [box_fun] = maxwell_box_smooth(center, box_size, varargin)


        %
        % Validate and parse inputs.
        %

    validateattributes(center, {'double'}, {'numel', 3, 'nonnan', 'finite'}, ...
                        mfilename, 'center');

    validateattributes(box_size, {'double'}, {'numel', 3, 'positive', 'finite'}, ...
                        mfilename, 'box_size');

    % Optional parameters.
    options = my_parse_options(struct('smooth_dist', 1), ...
                                varargin, mfilename);
    validateattributes(options.smooth_dist, {'numeric'}, ...
        {'positive', 'scalar', 'finite'}, mfilename, 'smooth_dist');


        %
        % Create function handle.
        %


        
    function [is_in, bounding_box] = f(x, y, z)
        bounding_box = {center - box_size/2 - options.smooth_dist, ...
                        center + box_size/2 + options.smooth_dist};

        s = options.smooth_dist;
        is_in = my_val_clamp(box_size(1)/2 - abs(x-center(1)), s) .* ...
                my_val_clamp(box_size(2)/2 - abs(y-center(2)), s) .* ...
                my_val_clamp(box_size(3)/2 - abs(z-center(3)), s);
    end

    box_fun = @f;
end


