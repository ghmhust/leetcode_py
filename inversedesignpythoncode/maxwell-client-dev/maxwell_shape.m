%% maxwell_shape
% Add shapes of constant material to the simulation domain.

%%% Syntax
%
% * |eps = maxwell_shape(grid, eps, eps_val, shape_fun)|
%   modifies the shape structure |eps| by inserting |eps_val| 
%   in the volume described by |shape_fun|.
%
% * |[eps, mu] = maxwell_shape(grid, [eps, mu], [eps_val mu_val], shape_fun)|
%   modifies both |eps| and |mu|.
%
% * |... = maxwell_shape(..., 'upsample_ratio', ratio)|
%   allows for an upsampling ratio of |ratio|, defaults to |ratio = 1|.
%
% * |... = maxwell_shape(..., 'f_avg', f_avg, 'f_rep', f_rep)|
%   allows for custom functions which determine the averaging function
%   for a grid point (|f_avg|) and how values of |eps| (and |mu|) 
%   are replaced (|f_rep|).

%%% Description
% |maxwell_shape| modifies the permittivity (and permeability, if desired)
% of the simulation domain.
% It does so by computing the approximate fill-fraction of a shape
% for every grid cell, and then inserting the appropriately scaled
% material value in that cell.
% As such, only shapes with constant epsilon and mu are supported.
%
% |maxwell_shape| is designed to be both simple and extensible at the same time.
% This is achieved in the following way:
%
% * |shape_fun| is a relatively simple function handle that 
%   only has to determine whether a point is inside or outside a shape.
%   This is to enable the user to easily define arbitrary shapes
%   without having to worry about grids, grid offsets, etc...
%   Please use |maxwell_box| and |maxwell_cyl| as templates for creating
%   your own |shape_fun|.
%
% * The |'f_avg'| and |'f_rep'| functions are provided so that the user
%   can customize how the average material parameter of the cell should
%   be calculated.
%   This level of extensibility is provided because the averaging of 
%   material parameters is still largely an unsolved problem.
%   By default, |'f_avg'| is a simple averaging function and
%   |'f_rep'| weights the new material value against the existing one
%   according to a cell's fill-fraction.
%   Since these functions are somewhat involved, the advanced user is
%   directed to the source code for additional information.
%
% As a last note, one of the disadvantages to this approach is that the removal
% of shapes is imperfect, meaning that writing the same shape twice,
% with different material values, will not be exactly equivalent
% to simply writing it once with the second material value only.
% There will be slight differences at the edge of the shape.
% 

%%% Source code
function [eps, mu] = maxwell_shape(grid, eps_mu, val, f, varargin)


        %
        % Validate and parse inputs.
        %

    my_validate_grid(grid, mfilename);

    [eps, mu] = my_split(eps_mu, grid.shape, {'eps', 'mu'}, mfilename);
    my_validate_field(eps, grid.shape, 'eps', mfilename);

    if isempty(mu)
        validateattributes(val, {'double'}, {'scalar', 'nonnan', 'finite'}, ...
                           mfilename, 'eps_val'); 
    else
        my_validate_field(mu, grid.shape, 'mu', mfilename);
        validateattributes(val, {'double'}, {'numel', 2, 'nonnan', 'finite'}, ...
                           mfilename, '[eps_val mu_val]'); 
    end

    validateattributes(f, {'function_handle'}, {}, mfilename, 'f');

    % Optional arguments.
    options = my_parse_options(struct(  'upsample_ratio', 1, ...
                                        'f_avg', @default_f_avg, ...
                                        'f_rep', @default_f_rep), ...
                                varargin, mfilename);
    validateattributes(options.upsample_ratio, {'numeric'}, ...
        {'positive', 'integer', 'scalar'}, mfilename, 'upsample_ratio');
    validateattributes(options.f_avg, {'function_handle'}, {}, ...
                        mfilename, 'f_avg');
    validateattributes(options.f_rep, {'function_handle'}, {}, ...
                        mfilename, 'f_rep');

    % Test if we can get a bounding box.
    % TODO: Check bounding box has non-zero (positive) volume.
    try 
        [~, bnd_box] = f(0, 0, 0);
        validateattributes(bnd_box, {'cell'}, {'numel', 2});
        validateattributes(bnd_box{1}, {'numeric'}, {'numel', 3});
        validateattributes(bnd_box{2}, {'numeric'}, {'numel', 3});
    catch
        bnd_box = {[-Inf -Inf -Inf], [Inf Inf Inf]};
    end

    % Test if we can give multiple points to f.
    try 
        out = f([0 1], [2 3], [4 5]);
        multipoint = true;
    catch
        multipoint = false;
    end


        %
        % Update material fields.
        %

    [eps_grid_pos, mu_grid_pos] = my_s2pos(grid); % Get grid positions.

    % Update components of epsilon.
    params = {bnd_box, f, options.upsample_ratio, multipoint, ...
                options.f_avg, options.f_rep};
    for k = 1 : 3
        eps{k} = my_update(eps_grid_pos{k}, eps{k}, k, val(1), params{:}); 
        if ~isempty(mu)
            mu{k} = my_update(mu_grid_pos{k}, mu{k}, k, val(2), params{:});
        end
    end


function [mat] = my_update(pos, mat, dir, val, box, f, ...
                            up_ratio, multipoint, f_avg, f_rep)
% Updates a single component of a field.

    % Determine box on which to evaluate f.
    for k = 1 : 3
        if box{1}(k) <= pos{k}(1)
            ind = 1;
        elseif box{1}(k) > pos{k}(end) % Not in the space.
            return
        else
            ind = max(find(pos{k} <= box{1}(k)));
        end
        s{1}(k) = ind;

        if box{2}(k) >= pos{k}(end)
            ind = length(pos{k});
        elseif box{2}(k) < pos{k}(1) % Not in the space.
            return
        else
            ind = min(find(pos{k} >= box{2}(k)));
        end
        s{2}(k) = ind;
        if s{1}(k) == s{2}(k) % 2D special case.
            s{2}(k) = s{2}(k) + 1;
        end
    end


    % Produce upsampled grid.
    for k = 1 : 3
        c{k} = zeros(up_ratio * (s{2}(k) - s{1}(k)), 1);
        cnt = 0;
        for l = s{1}(k) : s{2}(k)-1
            c{k}(cnt*up_ratio+[1:up_ratio]) = pos{k}(l) + ...
                (([0.5 : up_ratio] ./ up_ratio) .* (pos{k}(l+1)-pos{k}(l)));
            cnt = cnt + 1;
        end
    end


    % Obtain upsampled values.
    [x, y, z] = ndgrid(c{1}, c{2}, c{3});
    if multipoint
        inside_shape = f(x(:), y(:), z(:));
    else
        for k = 1 : numel(x)
            inside_shape(k) = f(x(k), y(k), z(k));
        end
    end
    inside_shape = reshape(inside_shape, size(x));

    if up_ratio == 1
        fill_fraction = inside_shape;
    else
        % Downsample results by averaging.
        for i = 1 : (s{2}(1) - s{1}(1))
            for j = 1 : (s{2}(2) - s{1}(2))
                for k = 1 : (s{2}(3) - s{1}(3))
                    fill_fraction(i, j, k) = f_avg(inside_shape(...
                                        (i-1)*up_ratio+[1:up_ratio], ...
                                        (j-1)*up_ratio+[1:up_ratio], ...
                                        (k-1)*up_ratio+[1:up_ratio]), dir);
                end
            end
        end
    end

    % Make changes to the material.
    m = mat(s{1}(1):s{2}(1)-1, s{1}(2):s{2}(2)-1, s{1}(3):s{2}(3)-1);
    mat(s{1}(1):s{2}(1)-1, s{1}(2):s{2}(2)-1, s{1}(3):s{2}(3)-1) = ...
        reshape(f_rep(val, fill_fraction(:), m(:), dir), size(m));


function [res] = default_f_avg(z, dir)
    res = mean(z(:));

function [res] = default_f_rep(val, ff, m, dir) 
    res = ff * val + (1-ff) .* m;
