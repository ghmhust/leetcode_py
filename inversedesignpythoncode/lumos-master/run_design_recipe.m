%% run_design_recipe
% Performs a design "recipe配方", which is a collection of calls to lumos().

%% Description
% Strings字符串 together multiple calls to lumos to get a design.

function run_design_recipe(problem_name, recipe_name, varargin)
     % run_design_recipe('wgc_te', 'r1', 'p0', p(:))
     % verify话，就要p0具体的值。如果是前面的优化部分的话，就不用输入了。p0在里面自有定义。这一项可以空。
     % (problem_name, recipe_name, varargin)
     % problem_name代表wgc_te.m等
     % recipe_name一般选择r1，后面会跳过verify步骤；如果要直接验证就选择verify
     % 参数输入范例run_design_recipe(problem_name, 'verify', 'p0', p(:));

    % Detect the 2D option.
    % strfind(S1,S2):寻找S2是否匹配S1，返回出现位置，没有出现返回空数组 
    % strrep(str1,str2,str3) ，进行字符串替换，区分大小写 ，把str1中所有的str2字串用str3来替换
    if ~isempty(strfind(problem_name, '2D'))
       % 二维
        flatten_option = true;
        
        exec_problem_name = strrep(problem_name, '2D', '2D');  
        % exec_problem_name = strrep(problem_name, '2D', ''); 
    else 
        % 三维
        flatten_option = false;
        exec_problem_name = problem_name;
    end

    % Function handle for generating the problem.
    gen_problem = eval(['@', exec_problem_name]);
    % eval（s）即 把字符串s的内容当作语句来执行
    % exec_problem_name指代wgc_te等
    
    % Set up the results directory for this recipe run.
    my_run_dir = [results_dir(), problem_name, '_', recipe_name, filesep];
    % 这里results_dir()代表lumos/results文件夹路径，results_dir()最后是目录分隔符\
    % 在results文件夹下，建立对应problem_name的文件夹，比如wgc_te_r1的文件夹，里面存有仿真过程中的场数据
 
     if isdir(my_run_dir)
           rmdir(my_run_dir,'s'); % Recursive remove递归删除.
     end
    mkdir(my_run_dir);
    % isdir判断输入（字符串）是否表示一个文件夹，A是一个文件夹，返回逻辑1（true），否则返回0（false）
    % mkdir（' fj '）, 表示在当前路径创建名为 fj 的文件夹
    % rmdir（'fl'）,表示删除当前路径下名为 fl 的文件夹，s表示移除指定文件夹及其文件夹内的所有内容
    % filesep用于返回当前平台的目录分隔符，Windows是反斜杠(\)，Linux是斜杠(/)
    
    function [p] = run_step(problem, params, step_name)
        my_step_name = [my_run_dir, 'step', step_name];
        fprintf('\nRunning step: %s\n', my_step_name);
        use_restart = false;
        % fprintf函数可以将数据按指定格式写入到指定的文本文件中
        % 数据的格式化输出：fprintf(fid,format,variables)，按指定的格式将变量的值输出到屏幕或指定文件
        % fid为文件句柄，若缺省，则输出到屏幕；format用来指定数据输出时采用的格式，%s输出字符串，\n换行
        % step_name是A的时候，对应global；是B的时候，对应local优化
        % params是对结构变量p值得优化方式，参数形式如下所示
  
  % p = run_step(problem, {'global', 'density', p, ...
  %                             [options.num_iters, 1e-3]}, 'A');
  
        % problem这个涉及到平面选择矩阵的产生S
%         % For no error override.
%         [z, p, vis] = lumos(my_step_name, problem, params{:}, ...
%                             'restart', use_restart);
%         return 
        while true
            try
                [z, p, vis] = lumos(my_step_name, problem, params{:}, ...
                                    'restart', use_restart);
               % params{:}，意味着参数长度可变，参数一=params{1}
                  break;
            catch exception
                fprintf(getReport(exception, 'extended'));
                fprintf('\n');
                use_restart = true;
                continue;
            end
        end
    end

    function [phi] = switch_to_phi(p)
        p = reshape(p, problem.design_area);
        phi = init_phi(p, [0 1], 1e-2, 0.5 * [1 -1]);
        phi = phi(:);
    end

    function [p] = my_phi2p(phi)
        phi = reshape(phi, problem.design_area);
        p = phi2p(phi, [0 1]);
    end

    function [phi] = reinit_phi(phi)
        p = my_phi2p(phi);
        phi = switch_to_phi(p);
    end
        
    % reinit_phi = @(phi) switch_to_phi(phi2p(reshape(phi, problem.design_area), [0 1]));

    % Log the diary.
    diary([my_run_dir, 'diary.txt']);
    diary on;
    % diary off;
    % diary off在最后一行，必须和diary on成对使用，
    % 如果只运行完diary on就退出调试，会导致rmdir指令无法删除文件，报错
    % diary的主要作用是将matlab工作过程中的全部屏幕文字和数据以文本的方式记录下来，成为一个工作日志
    switch recipe_name
        case 'r1'
            % Options structure
            options = struct(   'num_iters', 300, ...
                                'skip_A', false, ...
                                'skip_B', false, ...
                                'num_C_steps', 1, ...
                                'p0', 3/4);

            % Parse optional parameters.
            for k = 2 : 2 : length(varargin)
                options = setfield(options, varargin{k-1}, varargin{k});
            end

            p = options.p0;

            if flatten_option
                init_S_type = 'isolate'; % 隔离
            else
                init_S_type = 'alternate';
            end

            % Global optimization for 100 steps.
            if ~options.skip_A
                % Construct diagonalized problem对角化问题 for step A.
                problem = gen_problem({'flatten', flatten_option, ...
                                        'S_type', init_S_type});
              
                % 这里gen_problem输入即为wgc_te函数输入

                p = run_step(problem, {'global', 'density', p, ...
                                [options.num_iters, 1e-3]}, 'A');
            end

            % Local density optimization.
            if ~options.skip_B
                % Construct diagonalized problem for step B.
                problem = gen_problem({'flatten', flatten_option, ...
                                        'S_type', init_S_type, ...
                                        'size', 'small'});

                p = run_step(problem, ...
                            {'local', 'density', p, options.num_iters}, 'B');
            end

            % Non-diagonal problem for step C and on.
            problem = gen_problem({'flatten', flatten_option, ...
                                    'S_type', 'average', ...
                                    'size', 'small'});

            % Switch to level-set.
            phi = switch_to_phi(p);
            for i = 1 : options.num_C_steps 
                phi = run_step(problem, ...
                                {'local', 'level-set', reinit_phi(phi), ...
                                options.num_iters}, ['C', num2str(i)]);
            end

            % Run the verify step.
            p = my_phi2p(phi);
            run_design_recipe(problem_name, 'verify', 'p0', p(:));

        case 'verify'
            options = struct();
            % Parse optional parameters.
            for k = 2 : 2 : length(varargin)
                options = setfield(options, varargin{k-1}, varargin{k});
            end

            p = options.p0;

            % Count the number of field objectives in the original problem.
            ref_problem = gen_problem({'flatten', flatten_option, ...
                                    'S_type', 'average', ...
                                    'size', 'small'});

            for i = 1 : length(ref_problem.opt_prob)
                fobj = ref_problem.opt_prob(i).field_obj;
                num_fobj(i) = size(fobj.C, 2);
            end

            % Run the enlarged problem.
            problem = gen_problem({'flatten', flatten_option, ...
                                    'S_type', 'average', ...
                                    'size', 'large'});

            run_step(problem, {'local-no-move', 'density', p, 1}, 'V');

            % Get the results.
            fprintf('\nVerification results:\n');
            res = load([my_run_dir, 'stepV_state.mat']);
            for i = 1 : length(res.progress.out_power)
                fprintf('Mode %d: ', i);
                num_redundant = length(res.progress.out_power{i}) / num_fobj(i);
                for j = 1 : num_fobj(i)
                    data = res.progress.out_power{i}...
                                    ((j-1)*num_redundant+1:j*num_redundant);
                    fprintf('%1.4f (%1.4f)', mean(data), std(data));
                    if j ~= num_fobj(i)
                        fprintf(', ');
                    end
                end
                fprintf('\n');
            end
        otherwise
            error('Unkown recipe.');    
    end


    diary off;
end
