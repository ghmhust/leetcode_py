%% io
% Create a simple structure for cataloguing��Ϊ...��Ŀ����Ŀ¼���ࣩ inputs and outputs.

function [io_spec] = io(port, mode, power)
    io_spec = struct(   'port', port, ...
                        'mode', mode, ...
                        'power', power);
end
