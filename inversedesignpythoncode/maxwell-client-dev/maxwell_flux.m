%% maxwell_flux
% Electromagnetic power passing through a finite plane of the simulation.

%%% Syntax
%
% * |P = maxwell_flux(grid, [E H], plane_pos, plane_size)|
%   returns the integral of the Poynting vector over a plane
%   in the simulation domain.
%
% * |P = maxwell_flux(grid, [E H], [E1 H1])|
%   uses |[E1 H1]|, which is assumed to be the output from 
%   the |maxwell_wgmode| function (or something equivalent),
%   in order to define the plane and also filter |[E H]|.
%   This is useful for computing the amount of power in
%   only a single mode.
%   Note that this may fail if large evanescent fields are present.

%%% Description
% |maxwell_flux| calculates the integrated Poynting vector across
% a simulation plane. 
% It can alternatively be used to calculate the power in a single 
% waveguide mode, while filtering out the power in other modes.
%

%%% Source code
function [P] = maxwell_flux(grid, E_H, varargin)

        %
        % Validate and parse inputs.
        %

    
    my_validate_grid(grid, mfilename);

    [E, H] = my_split(E_H, grid.shape, {'E', 'H'}, mfilename);
    my_validate_field(E, grid.shape, 'E', mfilename);
    my_validate_field(H, grid.shape, 'H', mfilename);

    switch length(varargin)
        case 1
            E1_H1 = varargin{1};
            [E1, H1] = my_split(E1_H1, grid.shape, {'E1', 'H1'}, mfilename);
            my_validate_field(E1, grid.shape, 'E1', mfilename);
            my_validate_field(H1, grid.shape, 'H1', mfilename);

            [plane_pos, plane_size] = deal([], []);
        case 2
            plane_pos = varargin{1};
            plane_size = varargin{2};
            validateattributes(plane_pos, {'numeric'}, ...
                    {'nonnan', 'finite', 'numel', 3}, mfilename, 'plane_pos');
            validateattributes(plane_size, {'numeric'}, ...
                    {'nonnan', 'numel', 3}, mfilename, 'plane_size');

            [E1, H1] = deal([], []);
        otherwise
            error('Invalid number of input parameters.');
    end


        %
        % Determine prop_dir, plane_pos, and plane_size.
        %

    if ~isempty(E1) % Find prop_dir for filtered case.
        % First find the ranges of non-zero E1 values in every direction.
        for k = 1 : 3
            [ind{k}{1}, ind{k}{2}, ind{k}{3}] = ...
                    ind2sub(size(E1{k}), find(E1{k}(:) ~= 0)); 
        end
        for k = 1 : 3
            ind1{k} = [ind{1}{k}(:); ind{2}{k}(:); ind{3}{k}(:)];
            ind_range(k) = max(ind1{k}) - min(ind1{k});
        end
        ind_range(find(grid.shape == 1)) = inf; % Can't propagate out of 2D.

        % Deduce the propagation direction. 
        [~, prop_ind] = min(ind_range);

        % Obtain the position of the plane.
        pos = my_s2pos(grid); 
        perp_ind = mod(prop_ind, 3) + 1;
        prop_pos = pos{perp_ind}{prop_ind}(mode(ind1{prop_ind}));

        % Construct plane_pos and plane_size vectors.
        for k = 1 : 3
            grid_extent(k) = sum(real(grid.s_prim{k}));
            if isinf(grid_extent(k))
                grid_extent = realmax;
            end
        end
        grid_center = grid.origin + grid_extent;

        plane_pos = grid_center;
        plane_pos(prop_ind) = prop_pos;

        plane_size = 2 * grid_extent;
        plane_size(prop_ind) = +inf;
    end
     
    [p0, p1, prop_dir, prop_in_pos_dir] = my_find_plane(grid, plane_pos, plane_size);

    if ~isempty(E1) % For filtering case, take the whole plane.
        for k = 1 : 3
            if k ~= prop_dir
                p0(k) = 1;
                p1(k) = grid.shape(k);
            end
        end
    end


        %
        % Filter, if needed.
        %

    [vec, unvec] = my_vec(grid.shape);

    function [F] = my_project(F1, F2)
    % Project z1 onto (normalized) z2.
        F = unvec(vec(F2) * dot(vec(F1), vec(F2)) / norm(vec(F2))^2);
    end

    % Filter.
    if ~isempty(E1)
        E1{prop_dir} = 0 * E1{prop_dir};
        H1{prop_dir} = 0 * H1{prop_dir};
        E = my_project(E, E1);
        H = my_project(H, H1);
        % [E, H] = deal(E1, H1); % Debugging purposes only.
    end


        %
        % Cut out slice.
        %

    % Cut out the plane.
    for k = 1 : 3
        sp{k} = grid.s_prim{k}(p0(k):p1(k));
        sd{k} = grid.s_dual{k}(p0(k):p1(k));
        E{k} = E{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3));
        H{k} = H{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3));
    end


        %
        % Calculate the power.
        %

    P = my_calc_power(prop_dir, sp, sd, E, H);

    if ~prop_in_pos_dir
        P = -P;
    end

    if ~isempty(E1)
        for k = 1 : 3
            E1{k} = E1{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3));
            H1{k} = H1{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3));
        end
        if my_calc_power(prop_dir, sp, sd, E1, H1) < 0
            P = -P;
        end
    end
end

function [P] = my_calc_power(prop_dir, sp, sd, E, H)

    % Always need to ignore the field in the direction of propagation.
    a_dir = mod(prop_dir, 3) + 1;
    b_dir = mod(prop_dir+1, 3) + 1;
    [xa, ya, Ea, Ha] = deal(my_s2d(sd{a_dir}), ...
                            my_s2d(sp{b_dir}), ...
                            E{a_dir}, ...
                            H{b_dir});

    [xb, yb, Eb, Hb] = deal(my_s2d(sp{b_dir}), ...
                            my_s2d(sd{a_dir}), ...
                            E{b_dir}, ...
                            -H{a_dir});
    [xa, ya] = ndgrid(xa, ya);
    [xb, yb] = ndgrid(xb, yb);
    P = 0.5 * real( dot(xa(:) .* ya(:) .* Ea(:), Ha(:)) + ...
                    dot(xb(:) .* yb(:) .* Eb(:), Hb(:)));
end
