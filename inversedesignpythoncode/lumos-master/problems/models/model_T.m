%% model_T
% One port on the left and one port on the right.
% Made for beam resonators.

function [mode, vis_layer] = model_T(omega, in, out, wg_types, options)

%% Output parameters
% Fills in everything for mode structures, except for the in and out fields.
% At the same time, make the in and out fields easier to specify.

    % Basic dimensions.
    dims = [160 60 40];

    wg_dirs = {'+', '-'};
    wg_ypos = {dims(2)/2, dims(2)/2};

    for i = 1 : 2
        wg_options(i) = struct( 'type', wg_types{i}, ...
                                'dir', wg_dirs{i}, ...
                                'ypos', wg_ypos{i}); 
    end

    [mode, vis_layer] = metamodel_T(dims, omega, in, out, wg_options, options);
end
