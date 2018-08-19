function [wg] = build_io(ports, io, varargin)
     % tf = iscell(A),返回逻辑1（true）如果A是元胞数组，否则返回0（false）
    if ~iscell(io) % io不是元胞数组，将其变为元胞数组
        io = {io};
    end

    if ~isempty(varargin) % 判断一个数组是否是空的，如果A是一个空数组，返回逻辑1（true）；其他，则返回0（false）
        % 这里varargin不是空数组返回1
        make_redundant = true;
        redundancy_num = varargin{1};
    else
        make_redundant = false;
    end

    % Allow for a simple way to specify指定 a waveguide mode.
    for k = 1 : length(io)
        port = ports{io{k}.port};

        % NOTE: Distinguish区别 input and output by the number of elements
        % in the power field of io.
        dir = port.dir;
        switch length(io{k}.power)
            case 1 % Input
                shift = port.J_shift;
            case 2 % Output
                shift = port.E_shift;
                % Need to flip翻转 the dir in this case.
                if dir(2) == '+'
                    dir(2) = '-';
                else    
                    dir(2) = '+';
                end
            otherwise
                error('Invalid power.');
        end
        
        % Shift things appropriately.
        for l = 1 : 2
            pos{l} = port.pos{l} + shift;
        end

        % Build the waveguide spec.
        if ischar(io{k}.mode)
            mode_num = port.(io{k}.mode);
        elseif isnumeric(io{k}.mode)
            mode_num = io{k}.mode;
        else
            error('Invalid mode designation');
        end

        if ~make_redundant
            wg(k) = struct( 'type', port.type, ...
                            'power', io{k}.power, ...
                            'pos', {pos}, ...
                            'dir', dir, ...
                            'mode_num', mode_num);
        else
            % Used to get an average power output.
            for i = 1 : redundancy_num
                prop_dir = find(dir(1) == 'xyz');
                red_shift = +1 * (dir(2) == '-') + ...
                            -1 * (dir(2) == '+');

                red_pos = pos;
                red_pos{1}(prop_dir) = red_pos{1}(prop_dir) + (i-1) * red_shift;
                red_pos{2}(prop_dir) = red_pos{2}(prop_dir) + (i-1) * red_shift;
                wg((k-1)*redundancy_num + i) = ...
                    struct( 'type', port.type, ...
                            'power', io{k}.power, ...
                            'pos', {red_pos}, ...
                            'dir', dir, ...
                            'mode_num', port.(io{k}.mode));
            end
        end
    end
