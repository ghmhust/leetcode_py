%% example2_metalmode.m
% Solve for the resonant mode of a metallic resonator.

%%% Syntax
%
% * |[omega, E, H, grid, eps] = example2_metalmode()| 
%   runs the example in 3D and
%   returns the (complex) frequency of the eigenmode as |omega|.
%
% * |... = example2_metalmode('flatten', true)| 
%   runs the example in 2D.
%   This is very useful for quick tests.
%
% * |... = example2_metalmode('sim_only', true)| 
%   Only performs the initial simulation and 
%   does not perform the eigenmode solve.
%

%%% Source code
function [omega, E, H, grid, eps] = example2_metalmode(varargin)


        %
        % Parse inputs.
        %

    options = my_parse_options(struct(  'delta', 25, ...
                                        'flatten', false, ...
                                        'hires_delta', [3 3 6], ...
                                        'view_only', false, ...
                                        'sim_only', false), ...
                                varargin, mfilename);


        %
        % Create grid.
        %

    omega = 2*pi/940;
    x = -600 : options.delta : 600;
    y = -600 : options.delta : 600;
    z = -700 : options.delta : 500;
    if options.flatten
        x = 0;
    end

    hires_option = {[0 0 -100], [120 120 220], [options.hires_delta]};

    [grid, eps, ~, J] = maxwell_grid(omega, x, y, z, ...
                                        'hires_box', hires_option, ...
                                        'growth_rate', 1.1);

    % Structure constants.
    radius = 50;
    height = 200;

    gaas = 3.5^2;
    silver = -40-4i;
    air = 1;
    my_inf = 1e4;


        % 
        % Draw the GaAs and the silver shapes.
        %

    fprintf('Building shapes [%dx%dx%d grid]...\n', grid.shape);

    eps = maxwell_shape(grid, eps, silver, ...
                    maxwell_box([0 0 -height/2], [my_inf, my_inf, height]));
    eps = maxwell_shape(grid, eps, gaas, ...
                    maxwell_cyl([0 0 -height], radius, 2*height));
    eps = maxwell_shape(grid, eps, gaas, ...
                    maxwell_box([0 0 z(1)], [my_inf my_inf 2*(-z(1)-height)]));

    if options.view_only
        maxwell_view(grid, eps, [], 'y', [0 nan nan]);
        omega = nan;
        E = nan;
        H = nan;
        return
    end



        %
        % Point excitation in center for initial excitation.
        %

    c = round(grid.shape/2); % Center.
    J{2}(:, :, end-12) = 1;


        %
        % Solve initial simulation.
        %

    fprintf('Solving for initial field... ');
    [E, H] =  maxwell_solve(grid, eps, J); % Use this solution as an initial guess.

    maxwell_view(grid, eps, E, 'y', [nan nan 0]);
    
    if options.sim_only
        omega = grid.omega;
        return
    end


        % 
        % Find the eigenmode, using the previous solution field as initial guess.
        %

    [omega, E, H] =  maxwell_solve_eigenmode(grid, eps, E);
    % maxwell_view(grid, eps, E, 'y', [nan nan 0]);
    fprintf('\n');
end


