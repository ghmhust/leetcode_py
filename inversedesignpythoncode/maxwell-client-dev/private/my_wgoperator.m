%% Private wg_operator function.
function [A, get_wg_fields] = wg_operator(omega, s_prim, s_dual, epsilon, mu, ...
                                            prop_dir, prop_step, shape)
% Builds the operator (represented by matrix A), which defines the eigenmode 
% problem.
% Also provides the function get_wg_fields to obtain relevant parameters from
% the solution to the eigenmode problem.

    % Indices of the non-propagating directions.
    xdir = mod(prop_dir + 1 - 1, 3) + 1; 
    ydir = mod(prop_dir + 2 - 1, 3) + 1; 

    % Create matrices.
    xyz = 'xyz';
    Dx = deriv(xyz(xdir), shape);
    Dy = deriv(xyz(ydir), shape);

    % Stretched-coordinate parameters.
    [s_prim_x, s_prim_y] = ndgrid(s_prim{xdir}, s_prim{ydir});
    [s_dual_x, s_dual_y] = ndgrid(s_dual{xdir}, s_dual{ydir});

    % Build matrices.
    my_diag = @(z) spdiags(z(:), 0, numel(z), numel(z));
    Dfx = my_diag(1./s_dual_x) * Dx;
    Dbx = my_diag(1./s_prim_x) * (-Dx');
    Dfy = my_diag(1./s_dual_y) * Dy;
    Dby = my_diag(1./s_prim_y) * (-Dy');
    eps_yx = my_diag([epsilon{ydir}(:); epsilon{xdir}(:)]);
    inv_eps_z = my_diag(epsilon{prop_dir}.^-1);
    mu_xy = my_diag([mu{xdir}(:); mu{ydir}(:)]);
    inv_mu_z = my_diag(mu{prop_dir}.^-1);

    % Build operator.
    % Taken from Section 2 of Veronis, Fan, J. of Lightwave Tech., 
    % vol. 25, no. 9, Sept 2007.
    A = -omega^2 * eps_yx * mu_xy + ...
        eps_yx * [Dfy; -Dfx] * inv_eps_z * [-Dby, Dbx] - ...
        [Dbx; Dby] * inv_mu_z * [Dfx, Dfy] * mu_xy;

    % Build secondary operator to compute full h-field.
    v2h = @(beta, v)  [v; ((inv_mu_z * [Dfx, Dfy] * mu_xy * v) ./ (-i * beta))];

    % Build secondary operator to compute the error in the wave equation.
    my_zero = sparse(prod(shape), prod(shape));
    my_eye = speye(prod(shape));
    h_curl = @(beta)   [my_zero,        -i*beta*my_eye,  Dby; ...
                        i*beta*my_eye,  my_zero,       -Dbx; ...
                        -Dby,           Dbx,            my_zero];
    e_curl = @(beta)   [my_zero,        -i*beta*my_eye, Dfy; ...
                        i*beta*my_eye,  my_zero,        -Dfx; ...
                        -Dfy,           Dfx,            my_zero];
    eps = [epsilon{xdir}(:); epsilon{ydir}(:); epsilon{prop_dir}(:)];
    m = [mu{xdir}(:); mu{ydir}(:); mu{prop_dir}(:)];

    h_err = @(beta, h) norm(e_curl(beta) * ((h_curl(beta) * h) ./ eps) - ...
                        omega^2 * (m .* h)) / norm(h);
    e_err = @(beta, e) norm(h_curl(beta) * ((e_curl(beta) * e) ./ m) - ...
                        omega^2 * (eps .* e)) / norm(e);

    % Secondary operator to compute e-field.
    v2e = @(beta, v) (h_curl(beta) * v2h(beta, v)) ./ (i*omega*eps);

    % Secondary operator to compute j-field (excitation).
    n = prod(shape);
    v2j = @(v) [v(n+1:2*n); v(1:n); zeros(n, 1)];

    % Secondary operator to switch from a vector to the ordered field 
    % representation.
    rs = @(z) reshape(z, shape);
    to_field = @(z) {rs(z(1:n)), rs(z(n+1:2*n)), rs(z(2*n+1:3*n))};
    [~, rev_order] = sort([xdir, ydir, prop_dir]);
    reorder = @(f) {f{rev_order(1)}, f{rev_order(2)}, f{rev_order(3)}};
    vec2field = @(z) reorder(to_field(z));

    % Secondary operator that returns ordered fields.
    % The fields are also normalized so that the E- and H-fields are those
    % which give a Poynting vector of 1.
    % Also, subsequent solves should give the fields with the same phase
    % factor.
    % Lastly, J should be normalized to produce waveguide modes of power
    % 1 in both directions.
    function [E, H, J, E_err, H_err] = wg_fields(beta, v)
        % Obtain E- and H-fields (still in vector form).
        e = v2e(beta, v);
        h = v2h(beta, v);

        % Compute a normalization factor for power of 1.
        % Also ensure that subsequent solves yield waveguide modes of the same phase.

        % This calculates the total power of the mode as is,
        % and then uses the inverse root as the amplitude of the normalization factor.
        d_prim_x = my_s2d(s_prim_x);
        d_dual_x = my_s2d(s_dual_x);
        d_prim_y = my_s2d(s_prim_y);
        d_dual_y = my_s2d(s_dual_y);
        d_factor = [d_dual_x.*d_prim_y; d_prim_x.*d_dual_y];
        norm_amplitude = abs(0.5 * real(...
                dot(d_factor(1:n).*e(1:n), h(n+1:2*n)) + ...
                dot(d_factor(n+1:2*n).*e(n+1:2*n), -h(1:n))))^-0.5;
%                 dot(e(1:n), h(n+1:2*n)) + ...
%                 dot(e(n+1:2*n), -h(1:n))))^-0.5;

        % Use the element of the E-field with largest magnitude as a phase reference.
        % Actually, only use the first element...
        [~, ind] = max(abs(e));
        ind = 13;
        norm_angle = -angle(e(ind));

        % The combined normalization factor.
        norm_factor = norm_amplitude * exp(i * norm_angle);

        % Apply the normalization factor so that the fields produce 
        % Poynting vector of 1.
        e = norm_factor * e;
        h = norm_factor * h;
        v = norm_factor * v;

        % Fields in vector-field form.
        E = vec2field(e);
        H = vec2field(h);

        % Normalization factor for current excitation (some special sauce!).
        % Result from a maximum beta of pi/prop_step for discretized grids.
        nf_j = 1 + cos(real(beta)/2 * prop_step);
        v_factor = [d_factor(1:n)./d_prim_y; d_factor(n+1:2*n)./d_prim_x];
        v = v ./ prop_step;
        J = vec2field(nf_j * v2j(v));
        % d_factor(1)

        % Error in waveguide mode equation.
        % Note that this may increase because of the correction to the beta term,
        % especially for larger values of beta.
        E_err = e_err(beta, e);
        H_err = h_err(beta, h);
    end

    get_wg_fields = @wg_fields; % Function handle to get all the important info.

end % End of wg_operator private function.


%% Other private functions.
function [D] = deriv(dir, shape)
% Private function for creating derivative matrices.
% Note that we are making the forward derivative only.
% Also, we assume periodic boundary conditions.

    shift = (dir == 'xyz'); % Direction of shift.

    % Get the displaced spatial markers.
    my_disp = @(n, shift) mod([1:n] + shift - 1, n) + 1;
    [i, j, k] = ndgrid(my_disp(shape(1), shift(1)), ...
                        my_disp(shape(2), shift(2)), ...
                        my_disp(shape(3), shift(3)));

    % Translate spatial indices into matrix indices.
    N = prod(shape);
    i_ind = 1 : N;
    j_ind = i + (j-1) * shape(1) + (k-1) * shape(1) * shape(2);

    % Create the sparse matrix.
    D = sparse([i_ind(:); i_ind(:)], [i_ind(:), j_ind(:)], ...
                [-ones(N,1); ones(N,1)], N, N);
end % End of deriv private function.


