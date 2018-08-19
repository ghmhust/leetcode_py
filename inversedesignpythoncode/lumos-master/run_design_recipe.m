%% run_design_recipe
% Performs a design "recipe�䷽", which is a collection of calls to lumos().

%% Description
% Strings�ַ��� together multiple calls to lumos to get a design.

function run_design_recipe(problem_name, recipe_name, varargin)
     % run_design_recipe('wgc_te', 'r1', 'p0', p(:))
     % verify������Ҫp0�����ֵ�������ǰ����Ż����ֵĻ����Ͳ��������ˡ�p0���������ж��塣��һ����Կա�
     % (problem_name, recipe_name, varargin)
     % problem_name����wgc_te.m��
     % recipe_nameһ��ѡ��r1�����������verify���裻���Ҫֱ����֤��ѡ��verify
     % �������뷶��run_design_recipe(problem_name, 'verify', 'p0', p(:));

    % Detect the 2D option.
    % strfind(S1,S2):Ѱ��S2�Ƿ�ƥ��S1�����س���λ�ã�û�г��ַ��ؿ����� 
    % strrep(str1,str2,str3) �������ַ����滻�����ִ�Сд ����str1�����е�str2�ִ���str3���滻
    if ~isempty(strfind(problem_name, '2D'))
       % ��ά
        flatten_option = true;
        
        exec_problem_name = strrep(problem_name, '2D', '2D');  
        % exec_problem_name = strrep(problem_name, '2D', ''); 
    else 
        % ��ά
        flatten_option = false;
        exec_problem_name = problem_name;
    end

    % Function handle for generating the problem.
    gen_problem = eval(['@', exec_problem_name]);
    % eval��s���� ���ַ���s�����ݵ��������ִ��
    % exec_problem_nameָ��wgc_te��
    
    % Set up the results directory for this recipe run.
    my_run_dir = [results_dir(), problem_name, '_', recipe_name, filesep];
    % ����results_dir()����lumos/results�ļ���·����results_dir()�����Ŀ¼�ָ���\
    % ��results�ļ����£�������Ӧproblem_name���ļ��У�����wgc_te_r1���ļ��У�������з�������еĳ�����
 
     if isdir(my_run_dir)
           rmdir(my_run_dir,'s'); % Recursive remove�ݹ�ɾ��.
     end
    mkdir(my_run_dir);
    % isdir�ж����루�ַ������Ƿ��ʾһ���ļ��У�A��һ���ļ��У������߼�1��true�������򷵻�0��false��
    % mkdir��' fj '��, ��ʾ�ڵ�ǰ·��������Ϊ fj ���ļ���
    % rmdir��'fl'��,��ʾɾ����ǰ·������Ϊ fl ���ļ��У�s��ʾ�Ƴ�ָ���ļ��м����ļ����ڵ���������
    % filesep���ڷ��ص�ǰƽ̨��Ŀ¼�ָ�����Windows�Ƿ�б��(\)��Linux��б��(/)
    
    function [p] = run_step(problem, params, step_name)
        my_step_name = [my_run_dir, 'step', step_name];
        fprintf('\nRunning step: %s\n', my_step_name);
        use_restart = false;
        % fprintf�������Խ����ݰ�ָ����ʽд�뵽ָ�����ı��ļ���
        % ���ݵĸ�ʽ�������fprintf(fid,format,variables)����ָ���ĸ�ʽ��������ֵ�������Ļ��ָ���ļ�
        % fidΪ�ļ��������ȱʡ�����������Ļ��format����ָ���������ʱ���õĸ�ʽ��%s����ַ�����\n����
        % step_name��A��ʱ�򣬶�Ӧglobal����B��ʱ�򣬶�Ӧlocal�Ż�
        % params�ǶԽṹ����pֵ���Ż���ʽ��������ʽ������ʾ
  
  % p = run_step(problem, {'global', 'density', p, ...
  %                             [options.num_iters, 1e-3]}, 'A');
  
        % problem����漰��ƽ��ѡ�����Ĳ���S
%         % For no error override.
%         [z, p, vis] = lumos(my_step_name, problem, params{:}, ...
%                             'restart', use_restart);
%         return 
        while true
            try
                [z, p, vis] = lumos(my_step_name, problem, params{:}, ...
                                    'restart', use_restart);
               % params{:}����ζ�Ų������ȿɱ䣬����һ=params{1}
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
    % diary off�����һ�У������diary on�ɶ�ʹ�ã�
    % ���ֻ������diary on���˳����ԣ��ᵼ��rmdirָ���޷�ɾ���ļ�������
    % diary����Ҫ�����ǽ�matlab���������е�ȫ����Ļ���ֺ��������ı��ķ�ʽ��¼��������Ϊһ��������־
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
                init_S_type = 'isolate'; % ����
            else
                init_S_type = 'alternate';
            end

            % Global optimization for 100 steps.
            if ~options.skip_A
                % Construct diagonalized problem�Խǻ����� for step A.
                problem = gen_problem({'flatten', flatten_option, ...
                                        'S_type', init_S_type});
              
                % ����gen_problem���뼴Ϊwgc_te��������

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
