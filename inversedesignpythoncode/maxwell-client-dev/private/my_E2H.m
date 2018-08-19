function [H] = my_E2H(grid, mu, E)
% Calculate H from E.

    [vec, unvec] = my_vec(grid.shape);
    [~, A2] = my_functional_A(grid);

    H = unvec(A2(vec(E))); % Curl operator.

    % If mu does not equal to 1, make a simple correction.
    if ~isempty(mu)
        for k = 1 : 3
            H{k} = H{k} ./ (-1i * grid.omega * mu{k});
        end
    end

end
