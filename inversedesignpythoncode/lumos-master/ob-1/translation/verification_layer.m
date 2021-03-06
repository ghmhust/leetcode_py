%% verification
% Builds the structure and returns the modes and their performance.

%% Description
% 

function [modes] = verification_layer(opt_prob, z, varargin)

    N = length(opt_prob);

    %% Compute x, if needed.
    if isempty(varargin) % No x submitted.
        % Start all simulations.
        for i = 1 : N
            cb{i} = opt_prob(i).solve_A(z, opt_prob(i).phys_res.b(z));
            % 此处通过函数句柄调用solve_local.m/solve_maxwell.m进行波导FDFD仿真
        end
       
        % Wait for all simulations to complete.
        done = false * ones(N, 1);
        while ~all(done)
            for i = 1 : N
                [x{i}, done(i)] = cb{i}();
            end
        end

    else
        x = varargin{1};
    end

    %% Calculate relevant performance and visualization information.
    for i = 1 : N
        fobj = opt_prob(i).field_obj;

        % Raw magnitude of output power.
        outputs = fobj.C' * x{i};
        modes(i).raw_output_mag = [fobj.alpha, abs(outputs), fobj.beta]; 

        % Output angles.
        modes(i).output_angle = angle(outputs);

        % Output power in original units (of power).
        for j = 1 : size(fobj.C, 2)
            modes(i).output_power(j,:) = ...
                (modes(i).raw_output_mag(j,:)./norm(fobj.C(:,j))^2).^2;
        end

        % Calculate the physics residual.
        pr = opt_prob(i).phys_res;
        modes(i).phys_res_norm = norm(pr.A(z) * x{i} - pr.b(z)) / norm(pr.b(z));

        % Get epsilon.
        modes(i).epsilon = opt_prob(i).get_epsilon(z);
        
        % Get the E-field.
        modes(i).E = opt_prob(i).unvec(x{i});
    end

end
