function [A1, A2] = my_functional_A(grid)
% Here, A1 referes to curl on H-fields, and A2 refers to curl on E-fields.

    n = prod(grid.shape);
    [vec, unvec] = my_vec(grid.shape);

    % Inline function that takes (backward) derivatives of H.
    [spx, spy, spz] = ndgrid(grid.s_prim{1}, grid.s_prim{2}, grid.s_prim{3});
    function [f] = d1(f, dir)
        switch dir
            case 'x'
                f = (f - f([end, 1:end-1],:,:)) ./ spx;
            case 'y'
                f = (f - f(:,[end, 1:end-1],:)) ./ spy;
            case 'z'
                f = (f - f(:,:,[end, 1:end-1])) ./ spz;
            otherwise
                error('dir must be ''x'', ''y'', or ''z''.');
        end
    end

    % Inline function that takes (forward) derivatives of E.
    [sdx, sdy, sdz] = ndgrid(grid.s_dual{1}, grid.s_dual{2}, grid.s_dual{3});
    function [f] = d2(f, dir)
        switch dir
            case 'x'
                f = (f([2:end, 1],:,:) - f) ./ sdx;
            case 'y'
                f = (f(:,[2:end, 1],:) - f) ./ sdy;
            case 'z'
                f = (f(:,:,[2:end, 1]) - f) ./ sdz;
            otherwise
                error('dir must be ''x'', ''y'', or ''z''.');
        end
    end

    % Inline functions for the matrices.
    function [E] = A1_fun(H)
        H = unvec(H);
        E{1} = (d1(H{3}, 'y') - d1(H{2}, 'z'));
        E{2} = (d1(H{1}, 'z') - d1(H{3}, 'x'));
        E{3} = (d1(H{2}, 'x') - d1(H{1}, 'y'));
        E = vec(E);
    end

    function [H] = A2_fun(E)
        E = unvec(E);
        H{1} = (d2(E{3}, 'y') - d2(E{2}, 'z'));
        H{2} = (d2(E{1}, 'z') - d2(E{3}, 'x'));
        H{3} = (d2(E{2}, 'x') - d2(E{1}, 'y'));
        H = vec(H);
    end

    [A1, A2] = deal(@A1_fun, @A2_fun);
end
