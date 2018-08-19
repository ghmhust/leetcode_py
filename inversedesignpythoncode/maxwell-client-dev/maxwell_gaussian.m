%% maxwell_gaussian
% Excitation source for free-space Gaussian modes.

%%% Syntax
%
% * |J = maxwell_gaussian(grid, eps, plane_pos, plane_size, polarization, focal_length, beam_diameter)|
%   computes the current source for a Gaussian beam that is focused 
%   a distance |focal_length| away from |plane_pos|.
%   The Gaussian beam's polarization is determined via |polarization|,
%   which must be |'x'|, |'y'|, or |'z'|.
%   The Gaussian beam's diameter (full-width half-maximum) is determined 
%   by |beam_diameter|.
%   Similar to the |maxwell_wgmode| function, the excitation is provided 
%   at the finite plane located at |plane_pos|, 
%   which is of size |plane_size|.
%   One of the elements of |plane_size| must be either |+inf| or |-inf|
%   in order to denote the directionality of the desired waveguide mode.
%
% * |... = maxwell_gaussian(grid, [eps mu], ...)|
%   allows for |mu| not equal to 1.

%%% Description
% |maxwell_gaussian| allows the user to easily generate input Gaussian beams
% with chosen beam diameters and focal length.
% |maxwell_gaussian| is actually just a wrapper for the |maxwell_fsmode| function
% which is able to generate arbitrary free-space input excitation for the user.
% As such, it requires that the material parameters be uniform 
% across the excitation plane.
%

%%% Source code
function [J] = maxwell_gaussian(grid, eps_mu, plane_pos, plane_size, ...
                                polarization, focal_length, beam_diameter)


        %
        % Validate and parse inputs.
        %

    validateattributes(polarization, {'char'}, ...
                    {'scalar'}, mfilename, 'polarization');
    if ~any(polarization == 'xyz')
        error('Polarization must be either ''x'', ''y'', or ''z''.');
    end

    validateattributes(beam_diameter, {'numeric'}, ...
                    {'positive', 'finite', 'real'}, mfilename, 'beam_diameter');

    mode_fun = gaussian(plane_pos, beam_diameter, find(polarization == 'xyz'));
    J = maxwell_fsmode(grid, eps_mu, plane_pos, plane_size, mode_fun, ...
                        'focal_length', focal_length);
end

function [fun] = gaussian(center, fwhm, pol)
% Generate the mode function.
    sigma = fwhm / (2 * sqrt(2 * log(2)));

    function [E] = mode_fun(w, x, y, z)
    % The mode function (Gaussian in all directions of polarization pol).
        r = sqrt(   (x - center(1)).^2 + ...
                    (y - center(2)).^2 + ...
                    (z - center(3)).^2);
        E = (w == pol) .* ...
            1/(sigma*sqrt(2*pi)) * exp(-r.^2 / (2*sigma^2));
    end

    fun = @mode_fun;
    % fun = @(w, x, y, z) (w == pol) * ones(size(x));
end
        
