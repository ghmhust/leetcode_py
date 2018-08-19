%% example1_ring.m
% Excite the clock-wise mode of a ring resonator.

%%% Syntax
%
% * |[E, H, grid, eps] = example1_ring()| 
%   runs the example in 3D and
%   returns the E- and H-field (|E|, |H|),
%   as well as |grid| and |eps| which are useful for visualizing 
%   the result with |maxwell_view|.
%
% * |... = example1_ring('flatten', true)| 
%   runs the example in 2D.
%   This is very useful for quick tests.
%

%%% Source code
function [E, H, grid, eps] = example1_ring(varargin)


        %
        % Parse inputs.
        %

    options = my_parse_options(struct('flatten', false), ...
                                varargin, mfilename);


        %
        % Create grid.
        %

    % Make a grid for a wavelength of 1550 nm.
    if options.flatten
        [grid, eps] = maxwell_grid(2*pi/1.55, -3:0.05:3, -3:0.05:3, 0); % Use this for 2D.
        m = 2;
    else
        [grid, eps] = maxwell_grid(2*pi/1.55, -3:0.05:3, -3:0.05:3, -1.5:0.05:1.5);
        m = 1;
    end


        %
        % Setup the ring.
        %

    % Structure constants.
    height = 0.2;
    ring_radii = [1.4 1.0];
    si_eps = 13;
    air_eps = 1;

    % Draw ring.
    function [out] = custom_ring(x, y, z)
        r = sqrt(x.^2 + y.^2);
        out = (r < ring_radii(1) & r > ring_radii(2)) & (abs(z) < height/2);
    end

    eps = maxwell_shape(grid, eps, si_eps, @custom_ring); % Insert ring.

%     % Alternatively, use two cylinders
%     eps = maxwell_shape(grid, eps, si_eps, ...
%                         maxwell_cyl([0 0 0], ring_radii(1), height));
%     eps = maxwell_shape(grid, eps, air_eps, ...
%                         maxwell_cyl([0 0 0], ring_radii(2), height));


        %
        % Solve for initial excitation.
        %

    % Excitation for the fundamental mode (of the ring's waveguide).
    J = maxwell_wgmode(grid, eps, [0 mean(ring_radii) 0], [+inf 2 2], 'mode_number', m);

    fprintf('Initial excitation -- ');
    [E, H] =  maxwell_solve(grid, eps, J);

    maxwell_view(grid, eps, E, 'y', [nan nan 0], 'field_phase', inf); % Visualize the excited waveguide.
end
