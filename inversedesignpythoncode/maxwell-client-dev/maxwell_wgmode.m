%% maxwell_wgmode
% Excite a waveguide mode. Can also used to filter for a particular mode.

%%% Syntax
%
% * |[J, E, H, beta] = maxwell_wgmode(grid, eps, plane_pos, plane_size)|
%   computes the fundamental waveguide mode 
%   at the finite plane located at |plane_pos|, 
%   and of size |plane_size|.
%   One of the elements of |plane_size| must be either |+inf| or |-inf|
%   in order to denote the directionality of the desired waveguide mode.
%   The E- and H-fields of the mode are returned,
%   as well as the current excitation needed to excite the mode (|J|).
%   These field values actually consist of the entire vector-field.
%   Lastly, the wave-vector (|beta|) of the mode is also returned.
%   The fundamental mode is assumed.
%
% * |... = maxwell_wgmode(grid, [eps mu], ...)|
%   allows for |mu| not equal to 1.
%
% * |... = maxwell_wgmode(..., 'mode_number', m)|
%   returns the |m|-th order mode, where |m = 1| denotes the fundamental mode.
%   Defaults to 1.
%
% * |... = maxwell_wgmode(..., 'view', true)|
%   plots the fields of the waveguide mode.

%%% Description
% |maxwell_wgmode| computes the propagation mode for a nanophotonic waveguide structure
% including the wave-vector, E- and H-fields, as well as the current excitation
% needed for uni-directional excitation.
% The current excitation can be used for near-perfect excitation of a waveguide,
% while the mode fields can be used for filtering out powers at output waveguides.
%

%%% Source code
function [J, E, H, beta] = solve_wgmode(grid, eps_mu, plane_pos, plane_size, varargin)

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

    validateattributes(plane_pos, {'numeric'}, ...
                {'nonnan', 'finite', 'numel', 3}, mfilename, 'plane_pos');
    validateattributes(plane_size, {'numeric'}, ...
                {'nonnan', 'numel', 3}, mfilename, 'plane_size');
    if length(find(isinf(plane_size))) ~= 1
        error('plane_size must have exactly 1 element equal to either +inf or -inf.');
    end

    % Optional arguments.
    options = my_parse_options(struct(  'mode_number', 1, ...
                                        'view', false), ...
                                varargin, mfilename);
    validateattributes(options.mode_number, {'numeric'}, ...
                        {'positive', 'integer'}, mfilename, 'mode_number');
    validateattributes(options.view, {'logical'}, ...
                        {'binary'}, mfilename, 'view');


        %
        % Find plane (sub-grid) on which to solve for the eigenmode.
        %

    % Determine desired direction of propagation.
    [p0, p1, prop_dir, prop_in_pos_dir] = ...
                                    my_find_plane(grid, plane_pos, plane_size);

    % Cut out the bounded plane.
    sub_shape = p1 - p0 + 1;
    for k = 1 : 3
        sp{k} = grid.s_prim{k}(p0(k):p1(k));
        sd{k} = grid.s_dual{k}(p0(k):p1(k));
        e{k} = eps{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3));
        m{k} = mu{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3));
    end
    
        %
        % Build both real-only and full-complex versions of the operator
        % for the waveguide mode within the plane.
        %

    % Step size in propagation direction.
    prop_step = real(grid.s_dual{prop_dir}(p0(prop_dir))); 

    % Full complex operator.
    [A, get_wg_fields] = my_wgoperator(grid.omega, sp, sd, e, m, ...
                                        prop_dir, prop_step, sub_shape);

    for k = 1 : 3
        sp_r{k} = real(sp{k});
        sd_r{k} = real(sd{k});
        e_r{k} = real(e{k});
        m_r{k} = real(m{k});
    end

    % Real-only operator.
    A_r = my_wgoperator(real(grid.omega), sp_r, sd_r, e_r, m_r, ...
                                        prop_dir, prop_step, sub_shape);


        %
        % Solve for largest-magnitude eigenvalue of the real operator 
        % This is done in order to obtain the appropriate shift, 
        % from which we can calculate the most negative eigenvalues.
        %

    % Use the power iteration algorithm.
    n = size(A_r, 1);
    v = randn(n, 1);
    for k = 1 : 20 % 20 iterations should always be enough for an estimate.
        v = A_r * v;
    end
    ev_max = (v' * A_r * v) / norm(v)^2; % Rayleigh quotient.
    shift = abs(ev_max); % Shift works for both positive and negative ev_max.


        %
        % Solve for the desired eigenvector of the real operator
        % Taking the real operator, we a few of the most negative eigenmodes,
        % and then choose the one we are interested in.
        %

    % Shift the matrix and find the appropriate eigenmode.
    % Find a few extra modes just to be sure we found the correct one.
    [V, D] = eigs(A_r - shift * speye(n), options.mode_number + 2); 

    
    gamma = diag(D);
    [temp, ind] = sort(gamma); % Sort most negative first.
    v = V(:,ind(options.mode_number)); % Choose the appropriate eigenmode.


        %
        % Solve for the eigenvector of the full operator
        % We use the selected eigenvector from the real operator as an initial
        % guess.
        %

    % Perform Rayleigh quotient iteration to get the mode of the full operator.
    v_norqi = v;
    lambda = v' * A * v;
    rqi_done = false;
    for k = 1 : 40 
        err(k) = norm(A*v - lambda*v) / norm(lambda*v);
        if (err(k) < 1e-13)
            rqi_done = true;
            break
        end
        w = (A - lambda*speye(n)) \ v; 
        v = w / norm(w);
        lambda = v' * A * v;
    end
    if ~rqi_done 
        warning('Did not converge to mode, rqi error: %e.', err(end));
        v = v_norqi;
    end


        %
        % Calculate output parameters.
        %

    % Compute the wave-vector.
    beta = i * sqrt(lambda);
    beta = sign(real(beta)) * beta; % Force real part of beta to be positive.

    % Perform correction on beta to account for numerical dispersion.
    % Inspiration: Taflove's FDTD book, under Numerical Dispersion.
    % This correction term brings the error in emitted power to within 1 percent.
    % At the same time, additional error is introduced into the E_err and H_err terms.
    % This effect becomes more pronounced as beta increases.
    beta_corr = 2*sin(real(beta/2)) - real(beta);
    beta = beta + 0 * beta_corr; % Turn off correction.

    % Fields.
    [E_small, H_small, J_small, E_err, H_err] = get_wg_fields(beta, v);

    if prop_in_pos_dir
        coeff = -1;
    else
        coeff = +1;
    end

    % Make the components of the E and H fields match the propagation
    % direction.
    for k = 1 : 3
        if k ~= prop_dir
            H_small{k} = coeff * H_small{k};
        end
    end

    % Expand the fields to span the entire simulation space.
    for k = 1 : 3
        E{k} = zeros(grid.shape);
        H{k} = zeros(grid.shape);
        J{k} = zeros(grid.shape);

        E{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3)) = E_small{k};
        H{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3)) = H_small{k};
        J{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3)) = J_small{k};
    end

    %% Make the source uni-directional
    % This is done by creating an adjacent source which cancels the propagation
    % of the mode in one direction.

    dl = real(sp{prop_dir}); % Distance separating J and J_adj planes.


    % Shift indices for the propagation direction.
    ps0 = p0;
    ps1 = p1;
    ps0(prop_dir) = p0(prop_dir) + coeff;
    ps1(prop_dir) = p1(prop_dir) + coeff;

    % Form the adjacent J-field. 
    for k = 1 : 3  
        J{k}(ps0(1):ps1(1), ps0(2):ps1(2), ps0(3):ps1(3)) = ...
            -1 * J{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3)) * ...
            exp(-1i * beta * dl);
    end

    % Re-normalize so power flow is maintained at 1.
    for k = 1 : 3
        J{k} = J{k} ./ abs(1 - exp(coeff * 2i * beta * dl));
    end


        %
        % Plot fields, if desired.
        %

    if options.view
        f = {E_small{:}, H_small{:}};
        title_text = {'Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz'};
        for k = 1 : 6
            subplot(2, 3, k);
            my_plot(reshape(real(f{k}), sub_shape));
            title(title_text{k});
        end
        drawnow;
    end


function my_plot(x)
% Helps with plotting.
    if numel(find(size(x) ~= 1)) == 1 % Detect 1D data.
        plot([real(x(:)), imag(x(:))], '.-');
    else
        imagesc(squeeze(x).', (max(abs(x(:))) + eps) * [-1 1]);
        colorbar 
        axis equal tight;
        set(gca, 'YDir', 'normal');
    end
    colormap('jet');
