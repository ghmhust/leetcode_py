% TODO: Document, because this is tricky...
% Convention is that grid begins on primary grid points...
function [E_grid_pos, H_grid_pos] = my_s2pos(grid)
% Convert s-parameters to positions.
    origin_prim = grid.origin(:);
    origin_dual = grid.origin(:) + [real(grid.s_prim{1}(1))/2; ...
                                    real(grid.s_prim{2}(1))/2; ...
                                    real(grid.s_prim{3}(1))/2];

    for k = 1 : 3 % TRICKY!
        pos_prim{k} = origin_prim(k) + [0; cumsum(real(grid.s_dual{k}))];
        pos_dual{k} = origin_dual(k) + [0; cumsum(real(grid.s_prim{k}))];

        if any(isinf(pos_prim{k})) || any(isinf(pos_dual{k})) % 2D case.
            pos_prim{k} = grid.origin(k) * [1 1];
            pos_dual{k} = grid.origin(k) * [1 1];
        end
    end

    % Build up grid info.
    for k = 1 : 3
        for l = 1 : 3
            E_grid_pos{k}{l} = pos_prim{l};
            H_grid_pos{k}{l} = pos_dual{l};
        end
        E_grid_pos{k}{k} = pos_dual{k};
        H_grid_pos{k}{k} = pos_prim{k};
    end



