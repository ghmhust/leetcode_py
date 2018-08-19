function [A_diag] = my_A_diag(grid, eps_mu)

        %
        % Validate and parse input values.
        %

    my_validate_grid(grid, mfilename);

    [eps, mu] = my_split(eps_mu, grid.shape, {'eps', 'mu'}, mfilename);
    if isempty(mu)
        mu = my_default_field(grid.shape, 1); 
    end


        %
        % Compute the diagonal of the matrix for the wave equation in E.
        %

    
    [spx, spy, spz] = ndgrid(1./grid.s_prim{1}, 1./grid.s_prim{2}, ...
                                1./grid.s_prim{3});
    [sdx, sdy, sdz] = ndgrid(1./grid.s_dual{1}, 1./grid.s_dual{2}, ...
                                1./grid.s_dual{3});


    d{1} = (sdy + bs(sdy, 'y')) .* spy + (sdz + bs(sdz, 'z')) .* spz;
    d{2} = (sdz + bs(sdz, 'z')) .* spz + (sdx + bs(sdx, 'x')) .* spx;
    d{3} = (sdx + bs(sdx, 'x')) .* spx + (sdy + bs(sdy, 'y')) .* spy;

    for k = 1 : 3
        d{k} = d{k} - grid.omega^2 * eps{k};
    end

    A_diag = [d{1}(:); d{2}(:); d{3}(:)];


function [A] = bs(A, dir)
% Backward-shift of a matrix.
    switch dir
        case 'x'
            A = A([end, 1:end-1],:,:);
        case 'y'
            A = A(:,[end, 1:end-1],:);
        case 'z'
            A = A(:,:,[end, 1:end-1]);
        otherwise
            error('dir must be ''x'', ''y'', or ''z''.');
    end
