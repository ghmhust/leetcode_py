function [s_prim, s_dual] = my_pos2s(pos)

    for k = 1 : 3
        if length(pos{k}) == 1
            s_prim{k} = inf;
            s_dual{k} = inf;
        else
            s_prim{k} = my_prim(pos{k});
            s_dual{k} = diff(pos{k});
        end
    end


function [s] = my_prim(w)
% Private function to compute s_prim.
    w_avg = (w + [w(2:end); (w(end) + (w(2)-w(1)))]) ./ 2;
    s = diff(w_avg);
