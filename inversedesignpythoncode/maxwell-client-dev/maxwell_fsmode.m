%% maxwell_fsmode
% Excitation source for arbitrary free-space modes.

%%% Syntax
%
% * |J = maxwell_fsmode(grid, eps, plane_size, plane_pos, mode_fun)|
%   computes the current source needed to excite an arbitrary free-space mode.
%   The mode is determined via the user-supplied function handle |mode_fun|.
%   |mode_fun| must be of the form |E = mode_fun(w, x, y, z)| where
%   |E| is the E-field of the desired mode for a polarization of |w|
%   (which can be |1|, |2|, or |3| for x-, y-, or z-components of the E-field)
%   at the position |(x, y, z)|.
%   Similar to the |maxwell_wgmode| function, the excitation is provided 
%   at the finite plane located at |plane_pos|, 
%   which is of size |plane_size|.
%   One of the elements of |plane_size| must be either |+inf| or |-inf|
%   in order to denote the directionality of the desired waveguide mode.
%
% * |... = maxwell_fsmode(grid, [eps mu], ...)|
%   allows for |mu| not equal to 1.
%
% * |... = maxwell_fsmode(..., 'focal_length', flen)|
%   will attempt to produce the desired mode at a distance |flen|
%   away from the original excitation plane. 
%   If |flen| is non-zero, then the evanescent (non-propagating) components
%   are eliminated from the mode.
%   |flen| defaults to |0|.
%

%%% Description
% |maxwell_fsmode| allows the user to generate arbitrary free-space modes
% and even to produce them at a distance away from the excitation plane.
% It is able to do so by converting the mode into the plane-wave basis
% and back-propagating the mode by the appropriate distance.
% 
% Note that |maxwell_fsmode| requires that the material parameters be uniform 
% across the excitation plane. 
% Lastly, although the source is attempted to be made single-directional,
% this does not yet work quite as well as |maxwell_wgmode|, which is near perfect.
%

%%% Source code
function [J] = maxwell_fsmode(grid, eps_mu, plane_pos, plane_size, mode_fun, ...
                                varargin)


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
    options = my_parse_options(struct(  'focal_length', 0, ...
                                        'view', false), ...
                                varargin, mfilename);
    validateattributes(options.focal_length, {'numeric'}, ...
                        {'real', 'nonnan'}, mfilename, 'focal_length');
    validateattributes(options.view, {'logical'}, ...
                        {'binary'}, mfilename, 'view');


        %
        % Find plane (sub-grid) on which to put the free-space mode.
        %

    % Determine desired direction of propagation.
    [p0, p1, prop_dir, prop_in_pos_dir] = ...
                                    my_find_plane(grid, plane_pos, plane_size);

    % Get position information.
    pos = my_s2pos(grid);

    % Cut out the bounded plane.
    sub_shape = p1 - p0 + 1;
    for k = 1 : 3
        for l = 1 : 3
            pos{l}{k} = pos{l}{k}(p0(k):p1(k));
        end
        e{k} = eps{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3));
        m{k} = mu{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3));
    end

    % Get uniform values of epsilon and mu.
    eps = [e{1}(:); e{2}(:); e{3}(:)];
    mu =  [m{1}(:); m{2}(:); m{3}(:)];
    eps_val = mean(eps);
    mu_val = mean(mu);

    if any(eps ~= eps_val) || any(mu ~= mu_val)
        error('Material parameters must not vary over excitation plane.');
    end 
    
    % Step size in propagation direction.
    prop_step = real(grid.s_dual{prop_dir}(p0(prop_dir))); 

    % Helps us decide direction related signs.
    coeff = + 1 * (prop_in_pos_dir == true) - 1 * (prop_in_pos_dir == false);
    coeff = coeff * sign(options.focal_length);


        %
        % Get the mode shape.
        %

    % Test if we can give multiple points to f.
    try 
        out = mode_fun(1, [0 1], [2 3], [4 5]);
        multipoint = true;
    catch
        multipoint = false;
    end

    for k = 1 : 3
        [x, y, z] = ndgrid(pos{k}{1}, pos{k}{2}, pos{k}{3});
        if multipoint
            E{k} = mode_fun(k, x, y, z);
        else
            for cnt = 1 : numel(x)
                E{k}(cnt) = mode_fun(k, x(cnt), y(cnt), z(cnt));
            end
            E{k} = reshape(E{k}, size(x));
        end
    end

        
        %
        % Adjust to correct focal length.
        %

    omega_eff = real(grid.omega) * sqrt(real(eps_val) * real(mu_val));
    adj_prop_len = coeff * (abs(options.focal_length) + 0.5 * prop_step);
    E = my_propagate_beam(omega_eff, prop_dir, -adj_prop_len, E, pos);



        %
        % Form current source.
        %

    J = my_default_field(grid.shape, 0);

    if p0(prop_dir) ~= 1 && p0(prop_dir) ~= grid.shape(prop_dir)
        % Shifted positions for directionality
        ps0 = p0;
        ps1 = p1;
        ps0(prop_dir) = ps0(prop_dir) - coeff;
        ps1(prop_dir) = ps1(prop_dir) - coeff;

        for k = 1 : 3
            J{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3)) = E{k};
            if coeff ~= 0
                J{k}(ps0(1):ps1(1), ps0(2):ps1(2), ps0(3):ps1(3)) = -E{k} * ...
                        exp(-1i * prop_step * omega_eff);
            end 
            J{k} = J{k} ./ (1i * prop_step); % Scale with prop_step.
        end
    end
    

        %
        % Plot fields, if desired.
        %

    if options.view
        f = {E{:}};
        title_text = {'Ex', 'Ey', 'Ez'};
        for k = 1 : 3
            subplot(1, 3, k);
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
