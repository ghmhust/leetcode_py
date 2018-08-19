%% example4_focusedbeam
% Excite Gaussian or donut free-space modes.

%%% Syntax
%
% * |[E, H, grid, eps] = example4_focusedbeam('gaussian')| 
%   excites a Gaussian mode.
%
% * |... = example4_focusedbeam('donut')| 
%   excites a donut-mode (with radial H-field).
%
% * |... = example4_focusedbeam(..., 'flatten', true)| 
%   runs the example in 2D.
%   This is very useful for quick tests.
%
% * |... = example4_focusedbeam(..., 'flen', flen)| 
%   sets the focal length for the beam.
%   Defaults to |flen = 2|.
%

%%% Source code
function [E, H, grid, eps] = example4_focusedbeam(type, varargin)

        %
        % Parse inputs.
        %

    options = my_parse_options(struct(  'flatten', false, ...
                                        'flen', 2), ...
                                varargin, mfilename);


        %
        % Create grid.
        %

    % Make a grid for a wavelength of 1550 nm.
    omega = 2 * pi/1.55;
    x = -5 : 0.1 : 5;
    y = -5 : 0.1 : 5;
    z = -3 : 0.1 : 3;
    if options.flatten
        y = 0;
        % z = 0;
        mode_num = 1;
    else
        mode_num = 1;
    end


        %
        % Build simulation grid.
        %

    [grid, eps] = maxwell_grid(omega, x, y, z, ...
                            'hires_box', {[0 0 0], [1 .4 1], [.02 .02 .02]});

        %
        % Build excitation source.
        %

    switch type
        case 'gaussian'
            J = maxwell_gaussian(grid, eps, [0 0 2], [8 8 -inf], 'y', options.flen, 0.9);
            % J = maxwell_gaussian(grid, eps, [0 0 2], [8 8 -inf], 'y', options.flen, 0.9);
        case 'donut'
            mode_fun = zdonut([0 0 0], 0.8);
            J = maxwell_fsmode(grid, eps, [0 0 2], [8 8 -inf], mode_fun, 'focal_length', options.flen);
        otherwise
            error('Type must either ''gaussian'' or ''donut''.');
    end
   
    
        %
        % Solve simulation.
        %

    [E, H] =  maxwell_solve(grid, eps, J);
    maxwell_view(grid, eps, E, 'y', [nan 0 nan], 'field_phase', inf); % Visualize the excited waveguide.
end


function [fun] = zdonut(center, width)
% Generates the donut-mode excitation profile.

    function [E] = mode_fun(w, x, y, z)
    % Function for the donut-mode excitation.
        r = sqrt(   (x - center(1)).^2 + ...
                    (y - center(2)).^2) + 1e-10;
        E = (w == 1) .* (y./r) .* sin(pi*r/width/2) .* (r < 2*width) + ...
            (w == 2) .* (x./r) .* sin(pi*r/width/2) .* (r < 2*width);
    end

    fun = @mode_fun;
end
        
