function [p0, p1, prop_dir, prop_pos] = ...
                my_find_plane(grid, plane_pos, plane_size)

        %
        % Find plane (sub-grid) on which to solve for the eigenmode.
        %

    [eps_grid_pos, mu_grid_pos] = my_s2pos(grid); % Get grid positions.

    % Determine desired direction of propagation.
    prop_dir = find(isinf(plane_size));
    prop_pos = (plane_size(prop_dir) == +inf);

    % Determine the index of the plane (in propagation direction).
    avg_pos = mean([eps_grid_pos{prop_dir}{prop_dir}(:), ...
                    mu_grid_pos{prop_dir}{prop_dir}(:)], 2);
    [~, prop_ind] = min(abs(plane_pos(prop_dir) - avg_pos(1:end-1))); 

    % Determine plane for which to solve for the waveguide mode.
    pos = mu_grid_pos{3}; % Use the Hz point for reference.
    box = {plane_pos - plane_size/2, plane_pos + plane_size/2};
    for k = 1 : 3
        if k ~= prop_dir
            if length(pos{k}) == 2 % Special 2D case.
                p0(k) = 1;
                p1(k) = 1;

            else
                if box{1}(k) <= pos{k}(1)
                    ind = 1;
                else
                    ind = max(find(pos{k}(1:end-1) <= box{1}(k)));
                end
                p0(k) = ind;

                if box{2}(k) >= pos{k}(end-1)
                    ind = length(pos{k})-1;
                else
                    ind = min(find(pos{k}(1:end-1) >= box{2}(k)));
                end
                p1(k) = ind;
            end
        else
            p0(k) = prop_ind;
            p1(k) = prop_ind;
        end
    end
