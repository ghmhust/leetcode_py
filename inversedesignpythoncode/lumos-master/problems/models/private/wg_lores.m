%% wg_lores
% Provides nanophotonic waveguides for low-resolution designs.

%% Description
% Nanophotonic wavegudes are based on those from the litterature文献 which are designed
% to be single-mode at a wavelength around 1600 nm. 
%
% In this toolbox, the wavelength of 1600 nm corresponds 
% to the unit-less frequency of ~0.1571 (2*pi/40).

%
% For good performance, meaning minimal mode overap,
% the waveguide should reside in its own 40x40 box of SiO2.
%
% Assumes that the permittivity of silicon is 12.25,
% and that the cladding is SiO2 with a permittivity of 2.25.

function [waveguide, port] = wg_lores(epsilon, type, dir, len, pos)
    % [wg{i}, ports{i}] = wg_lores(epsilon, wg_options(i).type, ...
    %                            ['x', wg_options(i).dir], dims(1)-border, pos);

    % 三维情况：dims = [40 40 40];二维情况：dims = [40 40 1];
    % 介电常数数据结构：epsilon = {eps_lo*ones(dims), eps_lo*ones(dims), eps_lo*ones(dims)};
    % 分别代表xyz三个方向分量，每个分量在整个网格点上的值
    
    % size(a)表示矩阵每个维度的长度,size(a,n)表示矩阵a在第n个维度下的长度。
    % length(a)表示矩阵a的最大的长度，即max(size(a))
    % ndims(a)表示矩阵a的维数，即length(size(a))；
    % matlab认为向量也是二维矩阵，只不过其中一个维度的长为1.
    % squeeze函数用于删除矩阵中的维度值维1的那维，
    % 随机产生一个1x2x3的矩阵A，然后squeeze（A）将返回一个2x3的矩阵，将第一维却掉(因为第一位大小为1)
    eps_dims = size(epsilon{1});
    if ndims(epsilon{1}) == 2 % 二维
        flatten = true; 
        dims = [size(squeeze(epsilon{1})), 1];
    elseif eps_dims(3) == 1
        flatten = true;
        dims = [size(squeeze(epsilon{1})), 1];
    else 
        flatten = false; % 三维
        dims = size(epsilon{1}); % 网格点数目
    end
    
    % If a free-space wavelength is 40 grid points, 
    % and we want that to correspond to 1600 nm,
    % then the grid spacing is 40 nm.
    %grid_spacing = 10;
    grid_spacing = 40;
    
    % We assume a slab thickness of 250 nm.
    slab_thickness = 250 / grid_spacing;
    %slab_thickness = 250 / grid_spacing;
    
    wg_eps = 12.25;

    port.type = 'wgmode';
    port.dir = dir;

    switch type
        case 'single' % Contains a single TE and TM mode.
            width = 500 / grid_spacing;
            
            % Mode numbers.
            if flatten
                port.te0 = 2;
                port.tm0 = 1;
                % 二维
            else
                port.te0 = 1;
                port.tm0 = 2;
                % 三维
            end

        case 'double' % Contains up to first-order TE and TM modes.
            width = 650 / grid_spacing;

            % Mode numbers.
            if flatten
                port.te0 = 2;
                port.tm0 = 1;
                port.te1 = 4;
                port.tm1 = 3;
                % 二维
            else
                port.te0 = 1;
                port.tm0 = 2;
                port.te1 = 3;
                port.tm1 = 4;
                 % 三维
            end

        otherwise
            error('Unkown type.');
    end

    %port.pos = {pos - 19, pos + 20}; % Assumes 40x40 area...
    port.pos = {pos - 19, pos + 20}; % Assumes 40x40 area...
    switch dir(1)
        case 'x'
            wg_size = [len width];
            ind = 1;
        case 'y'
            ind = 2;
            wg_size = [width len];
        otherwise
            error('Unknown direction.');
    end
    port.pos{1}(ind) = pos(ind);
    port.pos{2}(ind) = pos(ind);

    % Make sure we haven't exceeded the limits of the space.
    for k = 1 : length(port.pos)
        port.pos{k} = (port.pos{k} < 1) * 1 + (port.pos{k} > dims) .* dims + ...
                    ((port.pos{k} >= 1) & (port.pos{k} <= dims)) .* port.pos{k};
    end


    % This is so that we are "centered" on the position given.
    % Centered means that E and J are in the correct positions and 
    % only use up +/- one plane from the position given.
    port.J_shift = zeros(1, 3);
    port.E_shift = zeros(1, 3);
    switch dir(2)
        case '+'
            port.J_shift(ind) = 0;
            port.E_shift(ind) = -1;
        case '-'
            % This looks weird, but it's because the waveguide mode solver
            % always adds one forward to compute a unidirectional mode.
            port.J_shift(ind) = -1;
            port.E_shift(ind) = +1;
        otherwise
            error('Unknown direction.');
    end

    center = size(epsilon{1})/2;
    waveguide = struct( 'type', 'rectangle', ...
                        'position', pos(1:2) - center(1:2), ...
                        'size', wg_size, ...
                        'permittivity', wg_eps);
