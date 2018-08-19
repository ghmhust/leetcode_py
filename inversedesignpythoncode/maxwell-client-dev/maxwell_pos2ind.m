%% maxwell_pos2ind
% Translate from position to array indices.

%%% Syntax
%
% * |[xi, yi, zi] = maxwell_pos2ind(grid, comp, pos)|
%   finds the array indices |[xi yi zi]| closest to |pos = [x y z]|
%   for the component |comp|.
%   |comp| must be either |'Ex'|, |'Ey'|, |'Ez'|, |'Hx'|, |'Hy'|, or |'Hz'|.
%
% * |... = maxwell_pos2ind(grid, comp, box)|
%   returns all the array indices within |box = {[x0 y0 z0], [x1 y1 z1]}|.
%
% * |[xi, yi, zi, ind] = maxwell_pos2ind(...)|
%   also returns |ind| which is the overall index corresponding to 
%   |[xi yi zi]|.
%

function [xi, yi, zi, ind] = maxwell_pos2ind(grid, comp, pos_or_box)

        %
        % Validate and parse input values.
        %

    my_validate_grid(grid, mfilename);

    validateattributes(comp, {'char'}, {'numel', 2}, 'comp', mfilename);
    complist = {'Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz'};
    if ~any(strcmp(comp, complist))
        error('Invalid parameter comp.');
    end

    if ~iscell(pos_or_box)
        box = {pos_or_box}; % Package pos as a cell array of 1.
    else
        box = pos_or_box;
    end
    validateattributes(box, {'cell'}, {'vector'}, ...
                        'box', mfilename);
    validateattributes(box{1}, {'numeric'}, ...
                        {'vector', 'numel', 3, 'nonnan', 'real'}, ...
                        'box{1}', mfilename);
    if length(box) == 2
        validateattributes(box{2}, {'numeric'}, ...
                            {'vector', 'numel', 3, 'nonnan', 'real'}, ...
                            'box{2}', mfilename);
    end


        %
        % Find the positions.
        %

    % Get positions.
    comp_ind = find(strcmp(comp, complist));
    [Epos, Hpos] = my_s2pos(grid);
    grid_pos = [Epos Hpos];
    pos = grid_pos{comp_ind};

    if length(box) == 1
        [~, xi] = min(abs(box{1}(1) - pos{1}));
        [~, yi] = min(abs(box{1}(2) - pos{2}));
        [~, zi] = min(abs(box{1}(3) - pos{3}));
        ind = sub2ind(grid.shape, xi, yi, zi);

    elseif length(box) == 2
        for k = 1 : 3
            in{k} = (pos{k} >= box{1}(k)) & (pos{k} <= box{2}(k));
        end
        [xi, yi, zi] = ndgrid(find(in{1}), find(in{2}), find(in{3}));
        ind = sub2ind(grid.shape, xi, yi, zi);
    end

    xi = xi(:);
    yi = yi(:);
    zi = zi(:);
    ind = ind(:);




