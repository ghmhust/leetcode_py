function [x] = my_insert_hires(x, center, width, delta, rate)
% Insert hi-resolution grid.

    if length(x) == 1 % Don't do anything for this case.
        return
    end

    % Determine index range to replace.
    ind0 = max([1; find(x <= center - width/2)]);
    ind1 = min([length(x); find(x >= center + width/2)]);

    % Form high-resolution region inside box.
    n = ceil((x(ind1) - x(ind0)) / delta);
    adjusted_delta = (x(ind1) - x(ind0)) / n;
    x_hires = x(ind0) + ([1:n-1] * adjusted_delta);

    % Form tapered wings.
    w1 = my_wing(x(ind1:end), adjusted_delta, rate);
    w0 = my_wing(-x(ind0:-1:1), adjusted_delta, rate);
    w0 = -w0(end:-1:1);

    % New positions.
    x = [w0(:); x_hires(:); w1(:)];


function [x] = my_wing(x, init_delta, rate)
% Taper off the high-resolution grid.
% Starts from x(1) and moves toward x(end).
% Assumes that we taper from high to low resolution.

    if length(x) == 1 % Don't need to add a wing if there's no room!
        return
    end

    if x(2)-x(1) <= init_delta*rate % No tapering needed.
        return
    end

    % Walk to find the end-point of the taper, 
    % as well as how many taper steps will be needed.
    cnt = 1;
    curr_x = x(1);
    while true
        curr_delta = init_delta * rate^cnt;
        curr_x = curr_x + curr_delta;

        ind = min(find(x >= curr_x));
        cnt = cnt + 1;

        if isempty(ind) % Hit end of structure.
            last_ind = length(x);
            wing_len = cnt + 1;
            break
        end

        if curr_delta >= x(ind) - x(ind-1) % Stop here, delta's "match".
            last_ind = ind - 1;
            wing_len = cnt + 1;
            break
        end

    end

%     if cnt == 2 % Special case where we shouldn't need tapering.
%         return
%     end

    % Form polynomial equation to solve for our adjusted rate.
    total_dist = x(last_ind) - x(1);
    c(1) = total_dist / init_delta;
    c(2) = - (init_delta + total_dist) / init_delta;
    c(wing_len+1) = 1;
    r = roots(c(end:-1:1));
    adjusted_rate = max(r(find(imag(r) == 0 & real(r) > 1)));

    if isempty(adjusted_rate) % Safety catch.
        adjusted_rate = rate;
    end

    % Form tapered wing.
    wing_x(1) = x(1);
    for k = 1 : wing_len-1
        wing_x(k+1) = wing_x(k) + init_delta * adjusted_rate^k;
    end

    % Replace values with taper.
    x = [wing_x(:); x(last_ind+1:end)];


