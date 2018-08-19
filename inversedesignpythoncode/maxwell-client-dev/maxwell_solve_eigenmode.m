%% maxwell_solve_eigenmode
% Eigenmode based on an initial guess for E.

%%% Syntax
%
% * |[omega, E, H] = maxwell_solve_eigenmode(grid, eps, E0)|
%   returns the frequency (|omega|), and fields (|E| and |H|)
%   of the eigenmode "nearest" to the field |E0|.
%   |maxwell_solve_eigenmode| works by calling 
%   an underlying |maxwell_solve| via the Rayleigh quotient iteration algorithm.
%
% * |... = maxwell_solve_eigenmode(grid, [eps mu], E0)|
%   does the same except for |mu ~= 1|.
%
% * |... = maxwell_solve_eigenmode(..., 'eig_max_iters', eig_n, 'eig_err_thresh', eig_err)|
%   sets the termination conditions for the eigenmode algorithm 
%   (Rayleigh quotient iteration).
%   Defaults to |eig_n = 20|, and |eig_err = 1e-9|.
%%
% * |... = maxwell_solve_eigenmode(..., 'vis_progress', vis_opt)|
%   controls the progress visualization for individual calls to |maxwell_solve|.
%   Defaults to |both|.
%
% * |... = maxwell_solve_eigenmode(..., 'max_iters', n, 'err_thresh', err)|
%   sets the termination conditions for the underlying calls to |maxwell_solve|.
%   Defaults to |n = 1e6| and |err = 1e-6|.

%%% Description
% |maxwell_solve_eigenmode| utilizes the Rayleigh quotient iteration algorithm
% to find an eigenmode of the system via repeated calls to |maxwell_solve|.
% To differentiate between the many eigenmodes which exist in the system,
% the user supplies an initial guess for the E-field of the eigenmode.
% In general, this guess is most easily produced by simulating
% the structure with a well-placed current excitation
% near the expected frequency of the electromagnetic mode.
%
% Finding electromagnetic modes of dispersive structures requires
% additional steps since material parameters which vary with frequency
% are not supported.
% A simple solution is to simply resolve the eigenmode with
% a structure tuned to the previously computed eigenmode frequency,
% and to iterate this way until the eigenmde frequency 
% no longer shifts significantly.
% Since PML's are actually dispersive, this can be done in the general case
% as well, by re-initializing the grid (using |'maxwell_grid'|) to the
% (real-part) of the previously computed eigenmode frequency.
%

%%% Source code
function [omega, E, H] = maxwell_solve_eigenmode(grid, eps_mu, E0, varargin) 

        %
        % Validate and parse inputs.
        %

    my_validate_grid(grid, mfilename);

    [eps, mu] = my_split(eps_mu, grid.shape, {'eps', 'mu'}, mfilename);
    if isempty(mu)
        mu = my_default_field(grid.shape, 1);
    end
    my_validate_field(eps, grid.shape, 'eps', mfilename);
    my_validate_field(mu, grid.shape, 'mu', mfilename);

    my_validate_field(E0, grid.shape, 'E0', mfilename);

    % Optional parameter-value pairs.
    options = my_parse_options(struct(  'eig_max_iters', 20, ...
                                        'eig_err_thresh', 1e-9, ...
                                        'vis_progress', 'both', ...
                                        'max_iters', 1e6, ...
                                        'err_thresh', 1e-6), ...
                                varargin, mfilename);


        %
        % Get ingredient matrices, vectors, and functions.
        %
        % The "F-field", which is defined as F = sqrt(epsilon) * E,
        % is used here because it transforms the wave equation from being
        % a generalized eigenvalue problem to a simple (canonical) one.
        % As such, most of the definitions in this section are in "F-space".
        %

    % Form the modified matrix and guess eigenvector. 
    [A0, v] = maxwell_axb(grid, [eps mu], E0, my_default_field(grid.shape, 0), ...
                            'functional', true);
    e = [eps{1}(:); eps{2}(:); eps{3}(:)];
    omega0 = grid.omega; % Constant omega for the inline function!!
    function [z] = A(z)
        z = A0(z) + omega0^2 * (e .* z);
    end
    v = v .* sqrt(e); % Transform to "F-space".

    % Compose function handles.
    mult_A = @(v) e.^-0.5 .* (A(e.^-0.5 .* v));
    % mult_A_dag = @(v) (e.^-0.5 .* (A.' * (e.^-0.5 .* conj(v)))).';
    sAinv_err = @(l, v, w) norm(mult_A(w) - l * w - v); % Useful for checking.

    % Helper functions.
    dims = grid.shape;
    n = prod(dims);
    unvec = @(z) {reshape(z(1:n), dims), reshape(z(n+1:2*n), dims), reshape(z(2*n+1:3*n), dims)};
    vec = @(z) [z{1}(:); z{2}(:); z{3}(:)]; 

    function [x] = solve_A_shifted(lambda, b)
    % Solves for the F-field.
        grid.omega = sqrt(lambda);
        J = unvec(sqrt(e) .* b ./ (-i * grid.omega));
        E = maxwell_solve(grid, [eps mu], J, ...
                            'max_iters', options.max_iters, ...
                            'err_thresh', options.err_thresh, ...
                            'vis_progress', options.vis_progress);
        x = sqrt(e) .* vec(E);
    end

    function my_vis(lambda, b, err)
        % Progress function.
        omega = sqrt(lambda);
        fprintf('wvlen: %1.3f, Q: %1.2e, omega: %1.1e + i%1.1e, err: %1.1e -- ', ...
                2*pi/real(omega), real(omega)/(2*imag(omega)/pi), real(omega), imag(omega), err);
    end


        %
        % Find the eigenmode.
        %

    [lambda, v] = my_solve_eigenmode(mult_A, @solve_A_shifted, @my_vis, v, ...
                                options.eig_max_iters, options.eig_err_thresh);
    fprintf('[eigenmode solve finished]\n');

    
        %
        % Back out the relevant parameters.
        %

    % Get back E and H.
    grid.omega = sqrt(lambda);
    omega = grid.omega;
    E = unvec(v ./ sqrt(e));
    H = my_E2H(grid, mu, E);

end
