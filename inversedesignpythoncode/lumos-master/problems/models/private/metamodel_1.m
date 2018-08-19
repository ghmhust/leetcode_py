%% metamodel_1
% Genericͨ�õ� in-plane nanophotonic model.

function [mode, vis_layer] = metamodel_1(dims, omega, in, out, ...
                                            wg_options, model_options)
%%
% wg_options�������������ͣ���ģ����ģ��������y��λ��
% model_optionsָ��model_options������size,flatten,s_type��Ϣ
%%
    %border = 60;% �߽�
    border = 15;% �߽�
    % ����border����verify����֤���׶Σ��ʵ�����ά�ȣ��÷�������Ľṹ�������п����ܷ����Ƿ���ƫ�ȷ����pml��ʹ����ȷ
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
    % dims(1:2)ָ��dims���1*3������ǰ����Ԫ��
    % boost���ӡ���ߣ�
    % large��ԭ�л��������ӷ������ߴ磬small����ԭ�з������ߴ粻��
   
    S_type = model_options.S_type;

    eps_lo = 2.25; % ��������
    eps_hi = 12.25; % ��
    z_center = dims(3)/2;
    z_thickness = 250 / 40; % �����
    %pml_thickness = [40 40 40]; % PML�߽�ߴ磬��λ40nm��Ҳ����������ʵ�ʳߴ��Ϊ400nm
    %z_thickness = 250 / 40; % �����
    pml_thickness = [10 10 10]; % PML�߽�ߴ磬��λ40nm��Ҳ����������ʵ�ʳߴ��Ϊ400nm
    % pml����dims�ڻ��ֵ�
    %reset_eps_val = eps_lo;
     reset_eps_val = eps_lo;
    % reset_eps_val�����������ʼ��糣��ֵ

    mu = {ones(dims), ones(dims), ones(dims)}; % ��ų���
    [s_prim, s_dual] = stretched_coordinates(omega, dims, pml_thickness);
    % the s-parameters which determine the spacing of the simulation grid
    % s_prim, s_dual��Ϊ1*3cell�������ݾ�Ϊ1*dims{1��2��3}����
    %% Construct structure
    epsilon = {eps_lo*ones(dims), eps_lo*ones(dims), eps_lo*ones(dims)};
    background = struct('type', 'rectangle', ...
                        'position', [0 0], ...
                        'size', [1e9 1e9], ...
                        'permittivity', eps_lo);
    % ע��background�����á������Ȼᱻ��д��ķ��������ǣ����Կ�ʼ���õĲ���x���򳤶Ⱥܳ���
    % metalmodel_1�����˳��Ⱥܳ��Ĳ�����Ȼ��add_planar���������background����д�븲�ǵ�һ���ֳ��ȵ������������
    % background��size��[1e9  1e9]�������ǵ����ַ����������������֮�󣬷�������С����
    % Ȼ��background�����������Ĳ��־�ֱ������
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
    % �������С
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
   
