%% example0_waveguide
% Excite the fundamental mode of a waveguide and measure output power.

%%% Syntax
%
% * |[E, H, grid, eps] = example0_waveguide()| 
%   runs the example in 3D and
%   returns the E- and H-field (|E|, |H|),
%   as well as |grid| and |eps| which are useful for visualizing 
%   the result with |maxwell_view|.
%
% * |... = example0_waveguide('flatten', true)| 
%   runs the example in 2D.
%   This is very useful for quick tests.
%

%%% Source code
function [E, H, grid, eps] = example0_waveguide(varargin)


        %
        % Parse inputs.
        %

    options = my_parse_options(struct(  'flatten', false), ...
                                varargin, mfilename);


        %
        % Create grid.
        %

    % Make a grid for a wavelength of 1550 nm.
    omega = 2 * pi/1.55;
    x = -2 : 0.025 : 2;
    y = -1 : 0.025 : 1;
    z = -1 : 0.025 : 1;
    if options.flatten
        z = 0;
        mode_num = 2;
    else
        mode_num = 1;
    end

    [grid, eps] = maxwell_grid(omega, x, y, z);


        %
        % Setup the waveguide.
        %

    % Structure constants.
    wg_height = 0.2;
    wg_width = 0.4;
    si_eps = 13;

    % Draw waveguide.
    eps = maxwell_shape(grid, eps, si_eps, ...
                        maxwell_box_smooth([0 0 0], [1e9 wg_width wg_height], ...
                                            'smooth_dist', 0.02));

        %
        % Solve for initial excitation.
        %

    [J, ~, ~, beta]  = maxwell_wgmode(grid, eps, [-1 0 0], [+inf 3 3], 'mode_number', mode_num);

    fprintf('Initial excitation -- ');
    [E, H] =  maxwell_solve(grid, eps, J);
    % �����ĳ��ֲ������в�����֧�ֵ�ģ���ĵ���
    
    % Visualize the excited waveguide.
    maxwell_view(grid, eps, E, 'y', [nan nan 0], 'field_phase', 0); 


        % 
        % Measure the output powers.
        %

    % Solve wgmode for filtering out the mode we care about.
    [~, E1, H1] = maxwell_wgmode(grid, eps, [0 0 0], [+inf 3 3], 'mode_number', mode_num);

    P0 = maxwell_flux(grid, [E H], [0 0 0], [+inf 100 100]);
    % ���㲨���������������ӡ͢ʸ���Ļ���
    P1 = maxwell_flux(grid, [E H], [E1 H1]);
    % �����ص����֣������������ģʽ�������ض�ģʽ����
    fprintf('Output powers at x = 0,\n');
    fprintf('Total power: %1.5f\n', P0);
    fprintf('Power in mode: %1.5f\n', P1);
