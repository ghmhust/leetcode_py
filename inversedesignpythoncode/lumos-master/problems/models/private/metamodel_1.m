%% metamodel_1
% Generic通用的 in-plane nanophotonic model.

function [mode, vis_layer] = metamodel_1(dims, omega, in, out, ...
                                            wg_options, model_options)
%%
% wg_options包括：波导类型（单模、多模）、方向、y轴位置
% model_options指代model_options，包括size,flatten,s_type信息
%%
    %border = 60;% 边界
    border = 15;% 边界
    % 增加border是在verify（验证）阶段，适当扩大维度，用仿真出来的结构，再运行看性能方面是否有偏差。确保，pml层使用正确
    if model_options.flatten % Make 2D.
        dims(3) = 1;
        pml_thickness(3) = 0;
    end

    if model_options.size == 'large'
         %size_boost = 80;
       size_boost = 20;
    elseif model_options.size == 'small'
        size_boost = 0;
    end
    dims(1:2) = dims(1:2) + 2 * size_boost;
    % dims(1:2)指代dims这个1*3向量的前两个元素
    % boost增加、提高，
    % large在原有基础上增加仿真区尺寸，small保持原有仿真区尺寸不变
   
    S_type = model_options.S_type;

    eps_lo = 2.25; % 二氧化硅
    eps_hi = 12.25; % 硅
    z_center = dims(3)/2;
    z_thickness = 250 / 40; % 硅层厚度
    %pml_thickness = [40 40 40]; % PML边界尺寸，单位40nm，也即三个方向实际尺寸均为400nm
    %z_thickness = 250 / 40; % 硅层厚度
    pml_thickness = [10 10 10]; % PML边界尺寸，单位40nm，也即三个方向实际尺寸均为400nm
    % pml是在dims内划分的
    %reset_eps_val = eps_lo;
     reset_eps_val = eps_lo;
    % reset_eps_val决定设计区初始介电常数值

    mu = {ones(dims), ones(dims), ones(dims)}; % 介磁常数
    [s_prim, s_dual] = stretched_coordinates(omega, dims, pml_thickness);
    % the s-parameters which determine the spacing of the simulation grid
    % s_prim, s_dual均为1*3cell，其内容均为1*dims{1或2或3}复数
    %% Construct structure
    epsilon = {eps_lo*ones(dims), eps_lo*ones(dims), eps_lo*ones(dims)};
    background = struct('type', 'rectangle', ...
                        'position', [0 0], ...
                        'size', [1e9 1e9], ...
                        'permittivity', eps_lo);
    % 注意background的作用。（长度会被新写入的仿真区覆盖，所以开始设置的波导x方向长度很长）
    % metalmodel_1设置了长度很长的波导，然后add_planar这个函数把background数据写入覆盖掉一部分长度的输入输出波导
    % background的size是[1e9  1e9]，它覆盖掉部分仿真区输入输出波导之后，仿真区大小不变
    % 然后background超出仿真区的部分就直接舍弃
    % Construct waveguides and modes.
    for i = 1 : length(wg_options)
        if wg_options(i).dir == '+'
            pos = [1+pml_thickness(1), wg_options(i).ypos+size_boost, z_center];
        elseif wg_options(i).dir == '-'
            pos = [dims(1)-pml_thickness(1)-1, wg_options(i).ypos+size_boost, z_center];
        else
            error('Unknown waveguide direction option');
        end
        [wg{i}, ports{i}] = wg_lores(epsilon, wg_options(i).type, ...
                                ['x', wg_options(i).dir], dims(1)-border, pos);
    end
    epsilon = add_planar(epsilon, z_center, z_thickness, {background, wg{:}});
 

    %% Build the selection matrix
    % Appropriate values of epsilon must be reset.
    design_pos = {border + [1 1] + size_boost, dims(1:2) - border - size_boost};
    design_area = design_pos{2} - design_pos{1} + 1;
    % 设计区大小
    [S, epsilon] = planar_selection_matrix(S_type, epsilon, ...
                                    design_pos, ...
                                    reset_eps_val, z_center, z_thickness);

      % Visualize the structure.
  for k = 1 : 3
      subplot(2, 3, k);
      imagesc(epsilon{k}(:,:)'); axis equal tight;
      % imagesc(epsilon{k}(:,:,20)'); axis equal tight;
      subplot(2, 3, k+3);
      imagesc(squeeze(epsilon{k}(:,40,:))'); axis equal tight;
      % imagesc(squeeze(epsilon{k}(:,40))'); axis equal tight;
  end
    %% Specify modes
    mode = struct(  'omega', omega, ...
                    'in', build_io(ports, in), ...
                    'out', build_io(ports, out), ...
                    's_prim', {s_prim}, ...
                    's_dual', {s_dual}, ...
                    'mu', {mu}, ...
                    'epsilon_const', {epsilon}, ...
                    'S', (eps_hi - eps_lo) * S, ...
                    'design_area', design_area);

    if model_options.size == 'large'
        % Add redundant out calculations.
        mode.out = build_io(ports, out, size_boost + 1);
    end


    %% Determine the visualization condition.
    if strcmp(in.mode(1:2), 'te')
        vis_component = 2; % Look at Ey.
    elseif strcmp(in.mode(1:2), 'tm')
        vis_component = 3; % Look at Ez.
    else
        error('Could not determine visualization component.');
    end

    vis_layer = struct( 'component', vis_component, ...
                        'slice_dir', 'z', ...
                        'slice_index', round(dims(3)/2));
   
