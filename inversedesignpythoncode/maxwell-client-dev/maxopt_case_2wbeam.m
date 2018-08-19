%% maxopt_case_2wbeam
% Sets up a frequency-doubling cavity optimization.

function [fun, x0] = maxopt_case_2wbeam(type, varargin)

        %
        % Parse inputs.
        %

    validateattributes(type, {'char'}, {'vector'}, 'type', mfilename);

    options = my_parse_options(struct(  'wvlen', [1550, 775], ...
                                        'eps_val', [3.2^2, 3.5^2], ...
                                        'delta', 40, ...
                                        'E', [], ...
                                        'flatten', false), ...
                                varargin, mfilename);


        %
        % Return recommended starting parameters.
        %

    
    N = 10;
    beam_height = 220;
    beam_width = 800;
    hole_pos = 350 * ([1:N] - 0.5);
    hole_xlen = 140 * ones(size(hole_pos));
    hole_ylen = (240 + 20*[1:N]) .* ones(size(hole_pos));

    hole_params = [hole_pos; hole_xlen; hole_ylen];
    x0 = [beam_height; beam_width; hole_params(:)];

    wvlen = options.wvlen;
    eps_val = options.eps_val; 


        %
        % Return appropriate function handle.
        %
        
    function [E, H, grid, eps] = get_fields(varargin)
        [~, ~, E, H, grid, eps] = solve_structure(varargin{:});
    end
        

    % Save previously used values.
    E_cache = options.E;

    function [fval, grad_f, omega, E, H, grid, eps] = cached_solve(varargin)
    % This is a cached solve which uses the most recent omega and E values
    % as initial guesses for the eigenmode solve.
        if isempty(E_cache)
            error('Need an initial E-field guess.');
        end

        [fval, grad_f, omega, E, H, grid, eps] = ...
            solve_eigenmodes(E_cache, varargin{:}); % Solve.

        % Update cached values.
        E_cache = E;
    end

    function [omega, E, H, grid, eps] = get_modes(varargin)
        [~, ~, omega, E, H, grid, eps] = cached_solve(varargin{:});
    end

    flt = options.flatten;
    switch type
        case 'get_fields'
            fun = @(x) get_fields(wvlen, eps_val, x, ...
                                    options.delta, flt, false);
        case 'fval'
            fun = @(x) solve_structure(wvlen, eps_val, x, ...
                                        options.delta, flt, false);
        case 'grad_f'
            fun = @(x) solve_structure(wvlen, eps_val, x, ...
                                        options.delta, flt, true);
        case 'get_fields_eig'
            fun = @(x) get_modes(wvlen, eps_val, x, ...
                                    options.delta, flt, false);
        case 'grad_f_eig'
            fun = @(x) cached_solve(wvlen, eps_val, x, ...
                                        options.delta, flt, true);
        otherwise
            error('Invalid type.');
    end
end

function [Fval, grad_F, omega, E, H, grid, eps] = ...
                solve_eigenmodes(E, wvlen, eps_val, varargin)
% Solve all eigenmodes.

    for k = 1 : length(wvlen)
        subplot(length(wvlen)+1, 1, k);
        [fval(k), grad_f{k}, omega{k}, E{k}, H{k}, grid{k}, eps{k}] = ...
                solve_one_eigenmode(E{k}, wvlen(k), eps_val(k), varargin{:});
    end
    subplot(length(wvlen)+1, 1, length(wvlen)+1);

    
%     % Optimize the lower-Q mode.
%     [Fval, ind] = max(fval);
%     grad_F = grad_f{ind};

%     % Favor the lower-Q mode.
%     Fval = 0;
%     grad_F = 0;
%     for i = 1 : length(grad_f)
%         if fval(i) < -1e4 % Q over 10k.
%             Fval = Fval + -1e4; % Cap fitness value at 10k for indiv. modes.
%             % Once Q exceeds 10k, no longer use that mode's gradient.
%         else % Normal addition.
%             Fval = Fval + fval(i);
%             grad_F = grad_F + grad_f{i};
%         end
%     end
%     ind = 0;

    % Another try.
    fmax = -1e4;
    Fval = 0;
    grad_F = 0;
    for i = 1 : length(fval)
        if fval(i) < fmax % Cap fitness to fmax.
            fval(i) = fmax;
        end
        fhat = abs((fval(i) - fmax) / fmax);
        Fval = Fval + 0.5*fhat^2;
        grad_F = grad_F + fhat * grad_f{i};
    end
    ind = 0;



    % Fval = sum(0.5 * (-1.0e-4 - fval).^2);

    % Pretty print.
    fprintf('fvals: ');
    for k = 1 : length(wvlen)
        if k == ind
            fprintf('[%e] ', fval(k));
        else
            fprintf('%e ', fval(k));
        end
    end
    fprintf('\n');
end


% Solve the resonance modes and return a gradient.
function [fval, grad_f, omega, E, H, grid, eps] = ...
                                solve_one_eigenmode(E, wvlen, eps_val, params, ...
                                                delta, flatten, calc_grad)


    [grid, eps, J] = make_structure(wvlen, eps_val, params, delta, flatten);


        %
        % Solve for the eigenmode.
        %

    [omega, E, H] = maxwell_solve_eigenmode(grid, eps, E, ...
                        'eig_max_iters', 10, 'vis_progress', 'text');
    maxwell_view(grid, eps, E, 'y', [nan nan 0], 'field_phase', 0); % Visualize.


        % 
        % Compose and evaluate fitness function (function to minimize).
        %

    function [fval, grad_w] = fitness(w)
    % Calculates figure of merit (fitness function) and its derivative.
        real_w = real(pi^2/wvlen);
        fval = -real_w / imag(w); % Negative Q-factor.
        grad_w = 1i * real_w / imag(w)^2;
    end
        
    [fval, grad_w] = fitness(omega);


        % 
        % Calculate structural gradient needed for gradient descent optimization.
        %

    if ~calc_grad % Skip if not needed.
        grad_f = nan;
        return
    end

    function [eps] = make_eps(params)
    % Function handle for creating the structure.
        [~, eps] = make_structure(wvlen, eps_val, params, delta, flatten);
    end

    function [lambda] = solver(eps)
    % Function that evaluates the fitness based on eps.
    % Only used for gradient checking.
        [omega_fit, E_fit, H_fit] = maxwell_solve_eigenmode(grid, eps, E, 'err_thresh', 1e-2);
        lambda = omega_fit^2;
    end

    % Calculate the structural gradient.
    grad_f = maxopt_freq_gradient(grid, E, omega, @fitness, params, @make_eps, ...
                'solver', @solver, ...
                'check_gradients', false);
end


function [Fval, grad_F, E, H, grid, eps] = ...
                solve_structure(wvlen, eps_val, varargin)
% Simulate all structures.
    % wvlen = wvlen(2);
    for k = 1 : length(wvlen)
        subplot(length(wvlen)+1, 1, k);
        [fval(k), grad_f{k}, E{k}, H{k}, grid{k}, eps{k}] = ...
                    solve_one_structure(wvlen(k), eps_val(k), varargin{:});
    end
    subplot(length(wvlen)+1, 1, length(wvlen)+1);

    % Find the worst one.
    [Fval, ind] = max(fval);
    grad_F = grad_f{ind};

    % Pretty print.
    fprintf('fvals: ');
    for k = 1 : length(wvlen)
        if k == ind
            fprintf('[%e] ', fval(k));
        else
            fprintf('%e ', fval(k));
        end
    end
    fprintf('\n');
end

function [fval, grad_f, E, H, grid, eps] = ...
                solve_one_structure(wvlen, eps_val, params, ...
                                    delta, flatten, calc_grad)
% Simulate a nanobeam structure.


    [grid, eps, J] = make_structure(wvlen, eps_val, params, delta, flatten);

        %
        % Use central point source as the excitation.
        %

    [x, y, z] = maxwell_pos2ind(grid, 'Ey', [0 0 0]);
    J{2}(x, y, z) = 1;
    J{2}(x-1, y, z) = 1;

        
        %
        % Solve.
        %

    [E, H] = maxwell_solve(grid, eps, J, 'vis_progress', 'text');
    maxwell_view(grid, eps, E, 'y', [nan nan 0], 'field_phase', 0); 


        %
        % Measure power reflected back to the center (figure of merit).
        %

    function [fval, grad_E] = fitness(E)
    % Calculates figure of merit and its derivative.
        % Figure of merit.
        E_meas = [E{2}(x, y, z); E{2}(x-1, y, z)];
        fval = mean((real(E_meas))); % This is the figure of merit.

        % Field gradient.
        grad_E = my_default_field(grid.shape, 0); 
        grad_E{2}(x, y, z) = 1/2;
        grad_E{2}(x-1, y, z) = 1/2;
    end
        
    [fval, grad_E] = fitness(E);


        % 
        % Calculate structural gradient.
        %

    if ~calc_grad % Skip if not needed.
        grad_f = nan;
        return
    end

    function [eps] = make_eps(params)
    % Function handle for creating the structure.
        [~, eps] = make_structure(wvlen, eps_val, params, delta, flatten);
    end

    % Calculate the structural gradient.
    grad_f = maxopt_field_gradient(grid, E, @fitness, params, @make_eps, ...
                'solver_fun', @(eps) maxwell_solve(grid, eps, J), ...
                'check_gradients', false, ...
                'solver_args', {'vis_progress', 'text'});
end


function [grid, eps, J] = make_structure(wvlen, eps_val, params, delta, flatten)

    height = params(1);
    width = params(2);
    hole_xpos = abs(params(3:3:end)); 
    hole_xlen = abs(params(4:3:end));
    hole_ylen = abs(params(5:3:end));


        %
        % Create grid.
        %

    % Make a grid for a wavelength of 1550 nm.
    omega = 2*pi / wvlen;
    x = -4200 : delta : 4200;
    y = -1200 : delta : 1200;
    z = -1000 : delta : 1000;

    if flatten
        z = 0;
    end

    [grid, eps, ~, J] = maxwell_grid(omega, x, y, z);


        %
        % Setup the structure.
        %

    % Structure constants.
    air_eps = 1;

    my_box = @(pos, siz) maxwell_box_smooth(pos, siz, 'smooth_dist', delta);
    
    % Draw the beam.
    eps = maxwell_shape(grid, eps, eps_val, ...
                        my_box([0 0 0], [1e9 width height]));

    % Draw rectangular holes.
    for k = 1 : length(hole_xpos)
        for l = [-1, 1]
            hole_pos = [l*hole_xpos(k) 0 0];
            hole_size = [hole_xlen(k) hole_ylen(k) 2*height];
            hole = my_box(hole_pos, hole_size);
            eps = maxwell_shape(grid, eps, air_eps, hole); 
        end
    end

    % Draw a central non-hole to make sure we don't cut out the center.
    eps = maxwell_shape(grid, eps, eps_val, my_box([0 0 0], [120 width height]-2*delta)); 
    % maxwell_view(grid, eps, [], 'y', [nan nan 0]);
end
