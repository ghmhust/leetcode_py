%% maxopt_example2_adjoint
% Form a cavity out of a square lattice using gradient optimization.


function [x, fval, f_vis] = maxopt_example2_adjoint(case_name, varargin)

    % case_name = 'squarepc';

        %
        % Parse inputs.
        %

    options = my_parse_options(struct(  'iters', 100, ...
                                        'case_args', {{}}, ...
                                        'max_delta', [], ...
                                        'init_step', [], ...
                                        'x0', [], ...
                                        'flatten', false), ...
                                varargin, mfilename);


        %
        % Set up the optimization problem.
        %

    flt = options.flatten;
    switch case_name
        case 'squarepc'
            [f, x0] = maxopt_case_squarepc('grad_f', 'flatten', flt);
            [f_vis] = maxopt_case_squarepc('get_fields', 'flatten', flt);
            init_step = 0.1;
            max_delta = 0.1;
        case '2wbeam'
            [f, x0] = maxopt_case_2wbeam('grad_f', 'flatten', flt, options.case_args{:});
            [f_vis] = maxopt_case_2wbeam('get_fields', 'flatten', flt, options.case_args{:});
            init_step = 1e2;
            max_delta = 10;
        case '2wbeam_eig'
           %  [f, x0] = maxopt_case_metalfocus('grad_f_eig', 'flatten', flt, options.case_args{:});
           % [f_vis] = maxopt_case_metalfocus('get_fields_eig', 'flatten', flt, options.case_args{:});
            [f, x0] = maxopt_case_2wbeam('grad_f_eig', 'flatten', flt, options.case_args{:});
            [f_vis] = maxopt_case_2wbeam('get_fields_eig', 'flatten', flt, options.case_args{:});
            init_step = 1e-1;
            max_delta = 10;
        case 'wdmgrating'
            [f, x0] = maxopt_case_wdmgrating('grad_f', 'flatten', flt);
            [f_vis] = maxopt_case_wdmgrating('get_fields', 'flatten', flt);
            init_step = 1e5;
            max_delta = 40;
        otherwise
            error('Invalid case_name.');
    end

    if ~isempty(options.init_step)
        init_step = options.init_step;
    end

    if ~isempty(options.max_delta)
        max_delta = options.max_delta;
    end

    if ~isempty(options.x0)
        x0 = options.x0;
    end

    % Visualization function for optimization progress.
    try
        mkdir(tempdir, case_name);
    end
    x_hist = [];
    function vis_progress(hist, step_size, x)
        fprintf('%d: %e [ss: %1.2e]\n', length(hist), hist(end), step_size); 
        plot(hist, '.-');
        xlabel('optimization iterations');
        ylabel('fval');
        title('structure optimization progress');
        saveas(gcf, [tempdir, case_name, filesep, ...
                    case_name, '_', sprintf('%04d', length(hist))], 'png');
        x_hist(:,length(hist)) = x(:);
        save([tempdir, case_name, filesep, 'x_hist.mat'], 'hist', 'x_hist');
    end

        
        %
        % Perform the optimization.
        %
        
    if options.iters > 0
        [x, fval, hist] = maxopt_gradient_descent(f, x0(:), ...
                                                    'init_step', init_step, ...
                                                    'max_delta', max_delta, ...
                                                    'max_iters', options.iters, ...
                                                    'vis_progress', @vis_progress);
    else
        x = x0;
        fval = nan;
    end

end
