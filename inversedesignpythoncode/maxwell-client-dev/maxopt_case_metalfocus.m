%% maxopt_case_metalfocus
% Used to optimize a metal focusing structure.

function [fun, x0] = maxopt_case_metalfocus(type, varargin)

        %
        % Parse inputs.
        %

    validateattributes(type, {'char'}, {'vector'}, 'type', mfilename);

    options = my_parse_options(struct(  'flatten', false), ...
                                varargin, mfilename);


        %
        % Return recommended starting parameters.
        %

    spacing0 = 10;
    y_roc = 400;
    z_roc = 400;
    z_depth = 250;
    x0 = [spacing0, y_roc, z_roc, z_depth];


        %
        % Return appropriate function handle.
        %
        
    function [E, H, grid, eps] = get_fields(varargin)
        [~, E, H, grid, eps] = solve_structure(varargin{:});
    end

    switch type
        case 'get_fields'
            fun = @(x) get_fields(x, options.flatten);
        case 'fval'
            fun = @(x) solve_structure(x, options.flatten);
        otherwise
            error('Invalid type.');
    end
end



function [fval, E, H, grid, eps] = solve_structure(params, flatten)
% Simulate the structure for the specific parameters..

        %
        % Get the structural parameters.
        %

    if isempty(params) % Used to look at the input beam.
        no_struct = true; 
    else
        no_struct = false;
        spacing = params(1); % Minimum spacing between the leads.
        y_roc = 2/params(2); % Radius of curvature in the y-direction.
        z_roc = 2/params(3); % Radius of curvature in the z-direction.
        z_depth = params(4); % How deep the silver penetrates the GaAs.
    end


        %
        % Create grid.
        %

    % Make a grid for a wavelength of 900 nm.
    delta = 30;
    wvlen = 900;
    omega = 2*pi/wvlen;

    my_size = 800;
    x = -my_size : delta : my_size;
    y = -my_size : delta : my_size;
    z = -my_size : delta : my_size;

    if flatten
        y = 0;
    end

    hibox = {[0 0 -z_depth/2], [50 50 50], [1 2 2]};
    [grid, eps, ~, J] = maxwell_grid(omega, x, y, z, 'hires_box', hibox);
%     [grid, eps, ~, J] = maxwell_grid(omega, x, y, z); 


        %
        % Setup the structure.
        %

    % Material permittivities.
    gaas = 3.5^2;
    silver = -40-4i;
    air = 1;

    % Draw GaAs slab.
    eps = maxwell_shape(grid, eps, gaas, ...
                        maxwell_box([0 0 min(z)], [inf inf 2*abs(min(z))]));

    % Draw hyperboloid leads.
    function [is_in] = hyperboloid_leads(x, y, z)
        is_in = (x.^2 - y_roc*y.^2 - z_roc*(z+z_depth/2).^2) >= spacing^2 & ...
                z <= 0 & z >= -z_depth;
        is_in = abs(x)-spacing - y_roc*y.^2 - z_roc*(z+z_depth/2).^2 >= 0 & ...
                z <= 0 & z >= -z_depth;
    end
    eps = maxwell_shape(grid, eps, silver, @hyperboloid_leads);

    subplot 121; maxwell_view(grid, eps, [], 'x', [nan 0 nan]);
    subplot 122; maxwell_view(grid, eps, [], 'x', [nan nan -z_depth/2]);


        %
        % Excite with a plane wave.
        %

    [~, ~, z] = maxwell_pos2ind(grid, 'Ex', [0 0 300]);
    J{1}(:,:,z) = 1;

    % Trick for 1-way plane-wave.
    J{1}(:,:,z) = 1 / grid.s_prim{3}(z);
    J_delta = grid.s_dual{3}(z);
    J{1}(:,:,z+1) = -exp(-1i * (2*pi/wvlen) * J_delta) / grid.s_prim{3}(z+1);
        

        %
        % Solve.
        %
    
    [E, H] = maxwell_solve(grid, eps, J);
    subplot 121; 
    maxwell_view(grid, eps, E, 'x', [nan 0 nan], 'field_phase', nan); 
    subplot 122; 
    maxwell_view(grid, eps, E, 'x', [nan nan -z_depth/2], 'field_phase', nan); 


        % 
        % Measure field strength near lead center.
        %

    [x, y, z] = maxwell_pos2ind(grid, 'Ex', [0 0 -z_depth/2]);
    fval = abs(E{1}(x, y, z));
end
