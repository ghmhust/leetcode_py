function [problem] = wgc_te2D(custom_model_options)
% 
%%
% custom_model_options采取如下所示注释里的形式(数据形式cell)
% problem = gen_problem({'flatten', flatten_option, ...
   %                                'S_type', 'average', ...
      %                              'size', 'small'});
%%
    % 功能：TE0转换为TE1
    % Choose the structure of the model (what waveguides to use where).
    model_structure = {'single', 'double'};
    % single和double分别代表500nm和650nm宽度的波导，也即单模和多模，详见wg_lores.m
    
    N = 1; % Number of modes.
%%
    % omega{1} = 2 * pi / 155;
    omega{1} = 2 * pi / 38.75;
    % 1550nm/40nm=38.75,其中40nm是网格精度，也即网格尺寸最小为40nm，其余长度依据网格精度进行划分
    in{1} = io(1, 'te0', 1);
    out{1} = {io(2, 'te1', [0.9 1]), io(2, 'te0', [0 0.01])};
    % 函数io(port, mode, power)，定义输入输出端口号、模式、功率范围
%% 
    vis_options.mode_sel = 1 : N;
    % 返回一维矩阵形式，其中N一般对应追踪模式的个数（不同波长、偏振。。。）
    
    % Build the problem.
    problem = get_problem(omega, in, out, vis_options, ...
                            @model_I, model_structure, custom_model_options);
    % 这里输入参量包括了  角频率、输入输出端口号、模式、功率范围、
    % 模式总数、对应模型、输入输出波导宽度类型、
    % 2D/3D仿真、S参数类型（average或alternate）、全空间size类型（small或large）
