%% stretched_coordinates
% Produces s-parameters needed to implement stretched-coordinate PML.

%% Description
% Produces the s-parameters needed to implement the stretched-coordinate 
% perfectly-matched layer (PML) boundary.
%
% Steven Johnson has a great reference on this.
%

function [s_prim, s_dual] = make_scpml(omega, origin, s_prim, s_dual, num_cells)

%% Input parameters
% * |omega| is the angular frequency of the simulation.
% * |t_pml| represents the depth, in grid points, of the pml in each direction.
%   It is also a 3-element vector. 
%   For no pml in a particular direction, set that element to 0.

%% Output parameters
% * |s_prim, s_dual| the primal and dual s-parameters for the grid.

%% Example
%
%   omega = 0.08;
%   dims = [80 40 20];
%   t_pml = [10 10 10];
%   [s_prim, s_dual] = stretched_coordinates(omega, dims, t_pml);

%% Source code

    % Position functions.
    pos = origin + [0; cumsum(s_dual(:))]; % Position of grid points.
    pos_prim = pos(1:end-1);
    pos_dual = (pos(1:end-1) + pos(2:end)) ./ 2;

    border = [pos(1) pos(end)];
    bnd = [pos(num_cells+1), pos(end-num_cells)];

    d = [bnd(1)-border(1), border(2)-bnd(2)];

    % Front PML.
    r = 1:num_cells;
    l_over_d = @(pos) (bnd(1) - pos) ./ (bnd(1) - border(1));
    s_prim(r) = stretch_s(omega, s_prim(r), l_over_d(pos_prim(r)), d(1));
    s_dual(r) = stretch_s(omega, s_dual(r), l_over_d(pos_dual(r)), d(1));

    % Back PML.
    r = length(s_prim) + [-num_cells+1:0];
    l_over_d = @(pos) (pos - bnd(2)) ./ (border(2) - bnd(2));
    s_prim(r) = stretch_s(omega, s_prim(r), l_over_d(pos_prim(r)), d(2));
    s_dual(r) = stretch_s(omega, s_dual(r), l_over_d(pos_dual(r)), d(2));
    return

function [s] = stretch_s(omega, s, l_over_d, d)
% We use the following formula for the imaginary part of s:
% -40 / (d * omega) * (l/d)^4
% where delta is the grid spacing, omega is the frequency,
% l is the distance inside the pml, and d is the width of the pml.
    s = s .* (1 - (40i * l_over_d.^4) ./ (d * real(omega)));
