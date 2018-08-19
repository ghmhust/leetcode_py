%% get_problem
% Generic���� function to create the design problem.

function problem = get_problem(omega, in, out, vis_options, ...
                                model, model_structure, custom_model_options)

    model_options = parse_model_options(custom_model_options);
   % ��ʼ��flatten��S_type��size��ֵ
   
    % Get the modes.��ȡģʽ��Ϣ
    for i = 1 : length(omega)
        [modes(i), vis_options(i).vis_layer] = ...
                                model(omega{i}, in{i}, out{i}, ...
                                        model_structure, model_options);
    end
    % ����model������Ӧ��model_I�Ⱥ�����������������������ʹ��ʱ��Ҫ����
    % length(a)��ʾ����a�����ĳ��ȣ���max(size(a))
    % vis_options(i).vis_layer��vis_options(i)����ṹ������
   
          % Visualize the structure.
  for k = 1 : 3
      subplot(2, 3, k);
      imagesc(modes.epsilon_const{k}(:,:)'); axis equal tight;  % 2D
      % imagesc(modes.epsilon_const{k}(:,:,20)'); axis equal tight;  % 3D
      subplot(2, 3, k+3);
      imagesc(squeeze(modes.epsilon_const{k}(:,40))'); axis equal tight;  % 2D
      % imagesc(squeeze(modes.epsilon_const{k}(:,40,:))'); axis equal tight;  % 3D
  end

    
    % Translate.
    if model_options.flatten
        solver = @solve_local; % 2D
    else
        solver = @solve_maxwell; % 3D
    end
    opt_prob = translation_layer(modes, solver);

    % Get design_area.
    design_area = modes(1).design_area;

    % Produce the final problem structure.
    problem = struct('opt_prob', opt_prob, ...
                    'vis_options', vis_options, ...
                    'design_area', design_area);
end
