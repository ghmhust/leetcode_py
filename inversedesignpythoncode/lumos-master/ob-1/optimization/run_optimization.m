%% run_optimization
% Obtain an optimized structure.

%% Description

function [z, p, state, progress] = ...
                run_optimization(opt_prob, g, p0, options, varargin)

    %% Set up logging.设置日志记录
    options = check_file_names(options);
    log_state = @(k, z, p, state, progress) save(options.state_file, ...
                        'k', 'z', 'p', 'state', 'progress');
    %log_state = @() save(options.state_file, ...
     %                   'k', 'z', 'p', 'state', 'progress');
    log_history = @(k, x, z, p) history_logger(options.history_file, ...
                                        {opt_prob.get_epsilon}, k, x, z, p);
    
    %% Initialize variables
    p = p0;
    z = g.m(p);
    state = [];
    progress = [];

    if ~isempty(varargin)
        state = varargin{1};
        progress = varargin{2};
    end

    k = options.starting_iter - 1;
    log_state(k, z, p, state, progress); 
    % save(options.state_file,'k', 'z', 'p', 'state', 'progress');
    % Log initial state.
    % 这一句报错，提示 错误使用load，无法读取文件‘E
    vp_state = [];

    termination_flag = false; % Delayed termination延迟终止 for local paradigm.

    prev_p = p0; % prev_p is used to revert back to previous structure.

    %% Run optimization
    for k = options.starting_iter : options.num_iters
        fprintf('%2d:', k);

        % Generate Q(z).
        if strcmp(options.paradigm, 'local')
            [P, q, state] = prodQ_local(z, opt_prob, state, ...
                                        options.paradigm_args{:});

            % Detect need to revert回复 back to previous structure.
            if ~state.prev_step_successful
                p = prev_p;
            end

            % Detect termination condition.
            if state.F == 0
                termination_flag = true;
            end
                                    
        elseif strcmp(options.paradigm, 'global')
            % Detect termination condition.
            if detect_global_termination(progress, options.err_thresh)
                fprintf(' [termination]\n');
                break
            end

            % Decide whether or not to reset dual variables u_i.
            if detect_u_reset(progress)
                fprintf(' [u_reset]');
                for i = 1 : length(state.u)
                    state.u{i} = 0 * state.u{i};
                end
                state.update_u = false;
            end

            % Produce Q(z).       
            [P, q, state] = prodQ_global(z, opt_prob, state, ...
                                        options.paradigm_args{:});
        else
            error('Invalid paradigm.');
        end

        % Early progress visualization.
        progress = options.vis_progress(k, state.x, z, p, progress);
        % options.vis_progress定义为函数句柄my_track_progress
        
        % Update the structure variable.
        if ~termination_flag % Only if not terminated.
            prev_p = p;
            [z, p] = update_structure(P, q, g, p, options.structure_args{:});  
        end

        % Log and visualize.
        progress = options.vis_progress(k, state.x, z, p, progress);
        log_state(k, z, p, state, progress);
        % save(options.state_file,'k', 'z', 'p', 'state', 'progress');
        log_history(k, state.x, z, p); % Do this last, in case hdf5 calls fail.

        fprintf('\n');

        if termination_flag% Termination for local paradigm.
            break
        end
    end
end % End run_optimization function.


%% Private functions
function [options] = check_file_names(options)
% Generate default names if needed for the state and history files.
    fname = strrep(datestr(clock), ' ', '@');
    if isempty(options.state_file)
        options.state_file = [tempdir, 'state_', fname, '.mat'];
        fprintf('Saving state data at: %s\n', options.state_file);
    end
    if isempty(options.history_file)
        options.history_file = [tempdir, 'hist_', fname, '.h5'];
        fprintf('Saving history data at: %s\n', options.history_file);
    end
end

function [reset_u] = detect_u_reset(progress)
% Decides whether or not the dual variables u should be reset.

    reset_u = false;
    if ~isstruct(progress)
        return
    end
    if length(progress.res_norm{1}) >= 2
        for i = 1 : length(progress.res_norm)
            prev_res_norm(i) = progress.res_norm{i}(end-1);
            curr_res_norm(i) = progress.res_norm{i}(end);
        end
        if max(curr_res_norm) > max(prev_res_norm)
            reset_u = true;
        end
    end
end

function [done] = detect_global_termination(progress, err_thresh)
% Decides whether or not the dual variables u should be reset.

    done = false;
    if ~isstruct(progress)
        return
    end

    for i = 1 : length(progress.res_norm)
        res_norm(i) = progress.res_norm{i}(end);
    end
    if max(res_norm) < err_thresh
        done = true;
    end
end

