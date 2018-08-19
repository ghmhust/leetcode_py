%% maxopt_gradient_descent
% Simple gradient descent optimization algorithm.

%%% Syntax
%
% * |maxopt_gradient_descent(fun, x0)|
%   attempts to minimize the function |fun| via gradient-descent
%   starting at the point |x0|.
%   |fun| is a function handle that returns both 
%   the value of the fitness function (figure of merit) and
%   the the gradient of said function with respect to |x|.
%   
% * |[x, fval, hist] = maxopt_gradient_descent(...)|
%   outputs |x|, the final 
%

%%% Description
% |maxopt_gradient_descent| is a simple function that performs a simple 
% gradient-descent optmization algorithm.
% |maxopt_gradient_descent| keeps a step even if it increases 
% (instead of decreasing) the figure of merit.
% However, in this case, the step-length will be decreased.
%

function [x_opt, f_opt, hist] = maxopt_gradient_descent(fun, x0, varargin)

        %
        % Validate and parse inputs.
        %

    validateattributes(fun, {'function_handle'}, {}, mfilename, 'fun');

    validateattributes(x0, {'numeric'}, ...
                        {'vector', 'real', 'finite', 'nonnan'}, ...
                        mfilename, 'x0');
        
    % Optional parameters.
    options = my_parse_options(struct(  'init_step', 1, ...
                                        'max_delta', 1, ...
                                        'step_shrink', 0.5, ...
                                        'step_big_shrink', 0.1, ...
                                        'step_grow', 1.1, ...
                                        'min_step', 1e-3, ...
                                        'max_iters', 100, ...
        'vis_progress', @(hist) fprintf('%d: %e\n', length(hist), hist(end))), ...
                                varargin, mfilename);

    simple_check = @(var, var_name) ...
        validateattributes(var, {'numeric'}, ...
                        {'positive', 'scalar', 'finite', 'nonnan'}, ...
                        mfilename, var_name);

    simple_check(options.init_step, 'init_step');
    simple_check(options.max_delta, 'max_delta');
    simple_check(options.step_shrink, 'step_shrink');
    simple_check(options.step_shrink, 'step_big_shrink');
    simple_check(options.step_grow, 'step_grow');
    simple_check(options.max_iters, 'max_iters');

    validateattributes(options.vis_progress, {'function_handle'}, {}, ...
                        mfilename, 'vis_progress');


        %
        % Perform minimization.
        %

    x = x0;
    [f, dx] = fun(x);
    step_size = options.init_step;
    f_opt = f;
    x_opt = x;
    for k = 1 : options.max_iters
        hist(k) = f;
        options.vis_progress(hist, step_size, x);

%         % fprintf('step size: %e\n', step_size);
%         max_step = step_size * max(abs(dx));
%         if max_step > options.max_delta
%             dx = dx ./ max_step;
%         end

%         dx_over = step_size * abs(dx) > options.max_delta;
%         dx = dx_over .* (dx ./ abs(dx)) * (options.max_delta * step_size) ...
%             + ~dx_over .* dx;

        x_step = step_size * dx;
        x_step = (abs(x_step) < options.max_delta) .* x_step + ...
                    (abs(x_step) >= options.max_delta) .* options.max_delta;
        x1 = x - x_step;
        [f1, dx1] = fun(x1);

        if f1 < f % Grow.
            step_size = step_size * options.step_grow;
        else % Shrink.
            step_size = step_size * options.step_shrink;
        end

        if (f1 < f/2) || (any([f f1]) > 0) % Accept as long as we still have half of original fitness.
            f = f1;
            x = x1;
            dx = dx1;
        else
            step_size = step_size * options.step_big_shrink;
        end

        if f1 < f_opt % Check if currently the best we got.
            f_opt = f1;
            x_opt = x1;
        end

        if step_size < options.min_step
            step_size = options.min_step; % Don't allow step_size below minimum.
        end
    end

        



