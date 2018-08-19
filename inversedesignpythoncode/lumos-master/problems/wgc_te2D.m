function [problem] = wgc_te2D(custom_model_options)
% 
%%
% custom_model_options��ȡ������ʾע�������ʽ(������ʽcell)
% problem = gen_problem({'flatten', flatten_option, ...
   %                                'S_type', 'average', ...
      %                              'size', 'small'});
%%
    % ���ܣ�TE0ת��ΪTE1
    % Choose the structure of the model (what waveguides to use where).
    model_structure = {'single', 'double'};
    % single��double�ֱ����500nm��650nm��ȵĲ�����Ҳ����ģ�Ͷ�ģ�����wg_lores.m
    
    N = 1; % Number of modes.
%%
    % omega{1} = 2 * pi / 155;
    omega{1} = 2 * pi / 38.75;
    % 1550nm/40nm=38.75,����40nm�����񾫶ȣ�Ҳ������ߴ���СΪ40nm�����೤���������񾫶Ƚ��л���
    in{1} = io(1, 'te0', 1);
    out{1} = {io(2, 'te1', [0.9 1]), io(2, 'te0', [0 0.01])};
    % ����io(port, mode, power)��������������˿ںš�ģʽ�����ʷ�Χ
%% 
    vis_options.mode_sel = 1 : N;
    % ����һά������ʽ������Nһ���Ӧ׷��ģʽ�ĸ�������ͬ������ƫ�񡣡�����
    
    % Build the problem.
    problem = get_problem(omega, in, out, vis_options, ...
                            @model_I, model_structure, custom_model_options);
    % �����������������  ��Ƶ�ʡ���������˿ںš�ģʽ�����ʷ�Χ��
    % ģʽ��������Ӧģ�͡������������������͡�
    % 2D/3D���桢S�������ͣ�average��alternate����ȫ�ռ�size���ͣ�small��large��
