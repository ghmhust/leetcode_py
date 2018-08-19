%% maxopt_case_wdmgrating
% Used to optimize a wavelength-splitting grating coupler.


function [fun, x0] = maxopt_case_wdmgrating(type, varargin)

        %
        % Parse inputs.
        %

    validateattributes(type, {'char'}, {'vector'}, 'type', mfilename);

    options = my_parse_options(struct(  'flatten', false), ...
                                varargin, mfilename);


        %
        % Return recommended starting parameters.
        %

    x0 = zeros(49, 2);
    if options.flatten
        x0 = zeros(7, 2);
    end
    wvlens = [1300 1500];
    N = length(wvlens);

%     [~, ~, E, H, grid, eps] = ...
%                 solve_structure(wvlens, ones(N, 1), x0, options.flatten, false);

        %
        % Calculate input powers.
        %

    if ~strcmp(type, 'get_fields')
        fprintf('Calculating input powers...\n');
        [~, ~, E, H, grid, eps] = ...
                    solve_structure(wvlens, ones(N, 1), [], options.flatten, false);
        for k = 1 : N
            P_in(k) = maxwell_flux(grid{k}, [E{k} H{k}], [0 0 0], [1e9 1e9 -inf]);
        end
    end


        %
        % Return appropriate function handle.
        %
        
    function [E, H, grid, eps] = get_fields(varargin)
        [~, ~, E, H, grid, eps] = solve_structure(varargin{:});
    end

    switch type
        case 'get_fields'
            fun = @(x) get_fields(wvlens, ones(1, N), x, options.flatten, false);
        case 'fval'
            fun = @(x) solve_structure(wvlens, P_in, x, options.flatten, false);
        case 'grad_f'
            fun = @(x) solve_structure(wvlens, P_in, x, options.flatten, true);
        otherwise
            error('Invalid type.');
    end
end


function [Fval, grad_F, E, H, grid, eps] = ...
            solve_structure(wvlens, P_in, params, flatten, calc_grad)

    N = length(wvlens);

    if flatten
        my_plane = [nan 0 nan];
    else
        my_plane = [nan nan 0];
    end

    % Initiate solves.
    for k = 1 : N
        subplot(1, N, k);
        [cb{k}, grid{k}, eps{k}, E_out{k}, J{k}] = ...
            start_structure_solve(wvlens(k), params, flatten);
    end

    % Wait for solves to finish.
    for k = 1 : N
        while ~cb{k}(); end
        [~, E{k}, H{k}] = cb{k}();

        % Visualize.
        subplot(1, N, k);
        maxwell_view(grid{k}, eps{k}, E{k}, 'y', my_plane);
    end

%     done = false * ones(N, 1);
%     while ~all(done)
%         for k = 1 : N
%             [done(k), E{k}, H{k}] = cb{k}();
%         end
%     end

    if isempty(params)
        [Fval, grad_F] = deal(nan);
        return
    end

    % General fitness function.
    [vec, unvec] = my_vec(grid{1}.shape);
    function [fval, grad_E] = fitness(E, E_ref, P_in)
        P_in = sqrt(abs(P_in) * norm(vec(E_ref))^2);
        overlap = vec(E_ref)' * vec(E);
        fval = (-1/P_in) * abs(overlap);
        grad_E = unvec((-1/P_in) * ((overlap)/abs(overlap)) * vec(E_ref));
    end
    
    % Compute fitness functions.
    fprintf('[');
    for k = 1 : N
        fitness_fun{k} = @(E) fitness(E, E_out{k}{k}, P_in(k));
        [fval(k), grad_E{k}] = fitness_fun{k}(E{k});
%         my_gradient_test(@(x) fitness_fun{k}(unvec(x)), ...
%                             vec(grad_E{k}), vec(E{k}), ...
%                             'real_with_imag', 'd');
        fprintf('%e ', fval(k));
    end
    [Fval, ind] = max(fval); % Find the worst performing wavelength.
    fprintf('\b]: %e\n', Fval);

    function [eps] = make_eps(params)
        [~, eps] = make_structure(wvlens(ind), params, flatten);
    end

    if calc_grad
        % Calculate gradient (if needed).
        grad_F = maxopt_field_gradient(grid{ind}, E{ind}, fitness_fun{ind}, ...
                    params, @make_eps, ...
                    'delta_p', 1e-3, ...
                    'solver_fun', ...
                            @(eps) maxwell_solve(grid{ind}, eps, J{ind}), ...
                    'check_gradients', false);
    else
        grad_F = nan;
    end

end
    
    

function [cb, grid, eps, E_out, J] = ...
                start_structure_solve(wvlen, params, flatten)
% Simulate the structure for the specific parameters.

    [grid, eps, J, wg_pos] = make_structure(wvlen, params, flatten);


        %
        % Initiate solve.
        %
    
    cb = maxwell_solve_async(grid, eps, J);


        %
        % Calculate output modes.
        %

    for k = 1 : 2
        if ~isempty(params)
            [J_wg{k}, E_out{k}, H_out{k}] = ...
                maxwell_wgmode(grid, eps, [2100 wg_pos(k) 0], [+inf 1e3 1e3]);
        else
            [J_wg{k}, E_out{k}, H_out{k}] = deal(nan);
        end
    end
    
%     subplot 121; 
%     maxwell_view(grid, eps, E, 'y', [nan -500 nan], 'field_phase', nan); 
%     subplot 122; 
%     maxwell_view(grid, eps, E, 'y', [nan nan 0], 'field_phase', nan); 
%     pause

end



function [grid, eps, J, wg_pos] = make_structure(wvlen, params, flatten)

        %
        % Get the structural parameters.
        %

    n = numel(params)/2;
    if n == 0 % Used to measure input power.
        no_struct = true;
    else
        no_struct = false;
        x_shift = params(1:n);
        y_shift = params(n+1:end);
        r_shift = zeros(n, 1); % No radius shift.
    end


        %
        % Create grid.
        %

    % Make a grid for a wavelength of 900 nm.
    delta = 40;
    omega = 2*pi/wvlen;

    my_size = 2500;
    x = -my_size : delta : my_size+500;
    y = -my_size : delta : my_size;
    z = -1000 : delta : 1000;

    if flatten
        y = -500;
    end

    [grid, eps] = maxwell_grid(omega, x, y, z); 


        %
        % Setup the structure.
        %

    % Material permittivities.
    si = 13;
    ox = 2.25;
    air = 1;

    % Structural constants.
    wg_pos = [-500 500];
    wg_width = 220;
    height = 440;
    radius = 150;
    a = 500;

    if ~no_struct
        % Draw the silicon slab. 
        eps = maxwell_shape(grid, eps, si, ...
                            maxwell_box([0 0 0], [4e3 4e3 height]));

        % Draw the silicon waveguides. 
        eps = maxwell_shape(grid, eps, si, ...
                maxwell_box([max(x) wg_pos(1) 0], [2*max(x) wg_width height]));
        eps = maxwell_shape(grid, eps, si, ...
                maxwell_box([max(x) wg_pos(2) 0], [2*max(x) wg_width height]));

        % Draw the holes.
        k = 1;
        j_range = -3 : 3;
        if flatten 
            j_range = -1;
        end
        for i = -3 : 3 
            for j = j_range
                pos = a * [i, j] + [x_shift(k), y_shift(k)];
                r = (radius + r_shift(k));
                if r > 0
                    my_cyl = maxwell_cyl_smooth([pos 0], r, 2*height, ...
                                                'smooth_dist', delta);
                    eps = maxwell_shape(grid, eps, air, my_cyl);
                end
                k = k + 1;
            end
        end

        % Draw oxide slab.
        eps = maxwell_shape(grid, eps, ox, ...
                    maxwell_box([0 0 min(z)-height/2+delta], ...
                                [inf inf 2*abs(min(z))]));
    end

%     subplot 121; maxwell_view(grid, eps, [], 'y', [nan wg_pos(1) nan]);
%     subplot 122; maxwell_view(grid, eps, [], 'y', [nan nan 0]);
%     pause

    if flatten
        my_plane = [nan 0 nan];
    else
        my_plane = [nan nan 0];
    end
    maxwell_view(grid, eps, [], 'y', my_plane);

        %
        % Excite with a Gaussian wave.
        %

    J = maxwell_gaussian(grid, eps, [0 0 500], [2*my_size 2*my_size -inf], ...
                        'y', 500, 2000);
end
