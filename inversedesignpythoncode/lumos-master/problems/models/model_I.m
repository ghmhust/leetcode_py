%% model_I
% One port on the left and one port on the right.
% ���Ҹ�һ���˿�

function [mode, vis_layer] = model_I(omega, in, out, wg_types, options)
%% 
% wg_typesָ��model_structure��single��ģ����500nm��ȣ�double��ģ����650nm
% optionsָ��model_options������size,flatten,s_type��Ϣ
%% Output parameters
% Fills in everything for mode structures, except for the in and out fields.
% At the same time, make the in and out fields easier to specify.

    % Basic dimensions.
    % dims = [360 280 160];
    dims = [90 70 40];
    % ����ȫ�ռ��С�����ϰ�����������͹裬�μ�����2016��Inverse Design of a Compact
    % Silicon-on-Insulator T-junction����ʾ��ͼ
    % �������½�Ϊ����ԭ�㣬
    % ��40nmΪ��λ��Ҳ��ʵ�ʳߴ�Ϊ3600��2800��1600nm
    wg_dirs = {'+', '-'};
    wg_ypos = {dims(2)/2, dims(2)/2};
   % x����Ⲩ��������y������ȷ���z������ȷ���
   % ������룬��Ϊ+���ұ������Ϊ-
   
    for i = 1 : 2
        wg_options(i) = struct( 'type', wg_types{i}, ...
                                'dir', wg_dirs{i}, ...
                                'ypos', wg_ypos{i}); 
    end
    % wg_typesָ��model_structure��single��ģ����500nm��ȣ�double��ģ����650nm
    % wg_options�������������͡�����y��λ��
    
    [mode, vis_layer] = metamodel_1(dims, omega, in, out, wg_options, options);
end
