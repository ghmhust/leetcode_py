% Allows for multiple simultaneous transfer over HTTP.
function my_stream_send (in, out)

    copier = MaxwellCopier; % Requires the Maxwell.jar library to be loaded.

    if length(in) ~= length(out)
        error('Unequal number of input and output streams.');
    end


    N = length(in);
    running = true * ones(N, 1);
    start_time = tic;
    status_time = start_time;
    prevlen = 0;

    while any(running)
        for k = 1 : N % Transfer some data.
            running(k) = copier.copy(in{k}, out{k});
        end

        if toc(status_time) > 0.3 || all(~running) % Periodically give updates.
            megabytes = copier.total_bytes_transferred / 1e6;
%             status_line = sprintf('%1.2f MB %s (%1.2f MB/s)', ...
%                 megabytes, action_name, megabytes/toc(start_time));
%             display_fun(status_line);
%             fprintf([repmat('\b', 1, prevlen), status_line]); % Write-over.
%             prevlen = length(status_line);
            status_time = tic;
        end
    end
    % fprintf('\n');
end

