%% model_I
% One port on the left and one port on the right.
% 左右各一个端口

function [mode, vis_layer] = model_I(omega, in, out, wg_types, options)
%% 
% wg_types指代model_structure，single单模波导500nm宽度，double多模波导650nm
% options指代model_options，包括size,flatten,s_type信息
%% Output parameters
% Fills in everything for mode structures, except for the in and out fields.
% At the same time, make the in and out fields easier to specify.

    % Basic dimensions.
    % dims = [360 280 160];
    dims = [90 70 40];
    % 仿真全空间大小，材料包括二氧化硅和硅，参见北邮2016年Inverse Design of a Compact
    % Silicon-on-Insulator T-junction论文示意图
    % 器件左下角为坐标原点，
    % 以40nm为单位，也即实际尺寸为3600、2800、1600nm
    wg_dirs = {'+', '-'};
    wg_ypos = {dims(2)/2, dims(2)/2};
   % x纵向光波传播方向、y波导宽度方向、z波导厚度方向
   % 左边输入，均为+，右边输出均为-
   
    for i = 1 : 2
        wg_options(i) = struct( 'type', wg_types{i}, ...
                                'dir', wg_dirs{i}, ...
                                'ypos', wg_ypos{i}); 
    end
    % wg_types指代model_structure，single单模波导500nm宽度，double多模波导650nm
    % wg_options包括：波导类型、方向、y轴位置
    
    [mode, vis_layer] = metamodel_1(dims, omega, in, out, wg_options, options);
end
