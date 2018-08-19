%% maxwell_solve_async
% Electromagnetic solve without waiting for completion.

%%% Syntax
% * |cb_fun = maxwell_solve_async(grid, eps, J)|
%   uploads a simulation and returns a callback function |cb_fun|
%   which is used to query progress查询进程 and return the results.
%   The callback function's syntax is |[done, E, H] = cb_fun()|
%   where |done| is a boolean variable布尔变量 that is set to true
%   only when |E| and |H| are the final solution fields.
%
% * |cb_fun = maxwell_solve_async(grid, [eps mu], J)|
%   is the same as above except that it allows |mu ~= 1|.
%
% The optional input parameters for |maxwell_solve| are also valid有效的 here.

%%% Description
% |maxwell_solve_async| is an asynchronous version 异步版本of |maxwell_solve|
% in that it exits as soon as一旦 the simulation has been uploaded to 
% the server, but before it is complete.
% For this reason, instead of returning the solution fields,
% a callback function is returned which is used to query查询 the
% progress of the simulation and to obtain the solution fields.
% The callback function can be used in this simple way (|maxwell_solve| does this):
%
%   while ~cb_fun(); end; % Wait until simulation finishes.
%   [~, E, H] = cb_fun(); % Get solution fields.
%
% The benefit of such a function is that it frees the local Matlab process
% to perform additional work.
%

%%% Source code
function [cb, vis_progress] = maxwell_solve_async(grid, eps_mu, J, varargin)

        %
        % Get current axis (for plotting) and start time (for timing).
        %

    my_axis = gca;


        %
        % Validate and parse inputs.
        %

    my_validate_grid(grid, mfilename);
    % mfilename是一个matlab内置的文件函数名，不用改。
    
    [eps, mu] = my_split(eps_mu, grid.shape, {'eps', 'mu'}, mfilename);
    if isempty(mu)
        mu = my_default_field(grid.shape, 1);
    end
    my_validate_field(eps, grid.shape, 'eps', mfilename);
    my_validate_field(mu, grid.shape, 'mu', mfilename);

    my_validate_field(J, grid.shape, 'J', mfilename);

    % Optional parameter-value pairs.
    options = my_parse_options(struct(  'E0', {my_default_field(grid.shape, 0)}, ...
                                        'vis_progress', 'both', ...
                                        'max_iters', 1e5, ...
                                        'err_thresh', 1e-6), ...
                                varargin, mfilename);
    my_validate_field(options.E0, grid.shape, 'E0', mfilename);
    validateattributes(options.vis_progress, {'char'}, ...
                {}, mfilename, 'vis_progress');
    validateattributes(options.max_iters, {'numeric'}, ...
        {'integer', 'positive', 'scalar', 'real'}, mfilename, 'max_iters');
    validateattributes(options.err_thresh, {'double'}, ...
        {'positive', 'scalar', '<', 1}, mfilename, 'err_thresh');


        %
        % Check if the simulation is 2D (can be solved locally).
        %
     
    no_print = false;
    function [varargout] = my_simple_callback(vis_progress, varargin)
        if strcmp(vis_progress, 'text') || strcmp(vis_progress, 'both') && ...
            ~no_print
            progress_text = '[finished] 2D problem solved locally';
            norm_p_text = [progress_text, ...
                    repmat(' ', 1, 60 - length(progress_text)), '\n'];
            fprintf(norm_p_text);
            no_print = true;
        end
        varargout = varargin;
    end
  
    if any(grid.shape == 1)
        % Compute E-field.
        [A, ~, b] = maxwell_axb(grid, [eps mu], options.E0, J);
        x = A \ b;
        N = prod(grid.shape);
        for k = 1 : 3
            E{k} = reshape(x((k-1)*N+1:k*N), grid.shape);
        end
        err = norm(A*x-b)/norm(b);
        
        % Compute H-field.
        H = my_E2H(grid, mu, E);

        % Return solution.
        vis_progress = options.vis_progress;
        cb = @() my_simple_callback(options.vis_progress, true, E, H, err);
        return
    end


        %
        % Upload simulation.
        %

    [server_url, name, vis_progress] = maxwell_upload(grid, eps, mu, J, ...
                                    options.E0, options.max_iters, ...
                                    options.err_thresh, options.vis_progress);


        %
        % Set up callback function.
        %

    % Persistent variables for the callback function.
    p_is_done = false;
    first_time = true;
    line_length = 60;
    p_E = [];
    p_H = [];
    p_err = [];
    p_state = [];
    no_print = false;
    
    start_time = tic;
    function [is_done, E, H, err] = maxwell_callback()
    % Queries server to inform user of the state of the simulation.
        if ~p_is_done % Not done, keep trying.
            [p_E, p_err, p_state, s] = maxwell_download(server_url, name);
            if strcmp(p_state, 'finished')
                p_is_done = true;
                p_H = my_E2H(grid, mu, p_E);
            end
        end

        % Update static variables.
        is_done = p_is_done;
        E = p_E;
        H = p_H;
        err = p_err;
        state = p_state;

        % Show the progress.
        last_print = false;
        if isempty(err) % Simulation not yet started.
            progress_text = sprintf('[%s] err: ----, iter: 0, seconds: %1.1f', ...
                                    state, toc(start_time));
        else 
            if is_done % Simulation complete.
                progress_text = sprintf('[%s] err: %e, iter: %d\n', ...
                                        state, err(end), length(err));
                last_print = true; % Make this the last line we print.

                
            else % Simulation in progress.
                progress_text = sprintf('[%s] err: %e, iter: %d, seconds: %1.1f', ...
                                        state, err(end), length(err), toc(start_time));
            end
        end

        if strcmp(vis_progress, 'text') | strcmp(vis_progress, 'both')
            % Normalized text progress output prints constant length of 60.
            norm_p_text = progress_text;
            if ~is_done
                norm_p_text = [norm_p_text, ...
                        repmat(' ', 1, line_length - length(progress_text))];
            end

            if ~first_time % If not first time, remove previous line.
                norm_p_text = [repmat('\b', 1, line_length), norm_p_text];
            end

            if ~no_print % Only print if this flag is false.
                fprintf(norm_p_text);
            end

            if last_print % No more printing!
                no_print = true;
            end

            first_time = false; % Denote that we've definitely printed once.
        end

        if strcmp(vis_progress, 'plot') | strcmp(vis_progress, 'both')
            % Plot the progress in log-scale.
            try
                axes(my_axis);
            catch 
                my_axis = gca;
                axes(my_axis);
            end

            if isempty(err)
                semilogy(1, 'bx');
            else 
                semilogy(err, 'b.-');
            end
            title(progress_text);
            xlabel('iterations');
            ylabel('error');

            % Add a dotted line showing the error threshold.
            hold on
            a = axis;
            semilogy(a(1:2), options.err_thresh * [1 1], 'k--');
            axis([a(1:2), options.err_thresh/10, a(4)]);
            hold off
        end

        % java.lang.System.gc(); % Request garbage collection for Java.
        pause(0.1);
    end

    cb = @maxwell_callback;
end


