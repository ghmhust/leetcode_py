%% maxwell_download
% Downloads result files from simulation server.

function [E, err, state, s] = maxwell_download(server_url, name)

   % modify_javapath();

    E = [];
    H = [];
    err = [];

    url = [server_url, name];
    endings = {'request', 'status', 'log', 'finished'};
    for k = 1 : length(endings)
        [s{k}, status(k)] = urlread([url, endings{k}]);
    end

    status_as_str = sprintf('%d%d%d%d', status); 
    switch status_as_str
        case '1000'
            state = 'queued';
        case '0010'
            state = 'loading';
        case '0110'
            state = 'solving';
        case '0111'
            state = 'finished';
        otherwise
            if status_as_str(4) == 1
                state = 'finished';
            else
                state = sprintf('unknown (%s)', status_as_str);
            end
    end

    if strcmp(state, 'solving') | strcmp (state, 'finished')
        % Get the error from the status file.
        err = str2num(s{2});
    end

    if strcmp(state, 'finished') 
        % Download all files.
        [quad, comp] = ndgrid('ri', 'xyz');
        for k = 1 : numel(quad)
            files{k} = [name, 'E_', comp(k), quad(k)];
            % copyfile('[server_url,files{k}]',server_url);
        end
        
        % my_download(files, tempdir, server_url);

        % Load º”‘ÿfiles.
        for k = 1 : 3
            E{k} = double(h5read([tempdir, files{2*k-1}], '/data')) + ...
                1i * double(h5read([tempdir, files{2*k}], '/data'));
            E{k} = permute(E{k}, [ndims(E{k}):-1:1]); % Convert to column-major.
        end

        % Delete files.
%       for k = 1 : numel(files)
 %           delete([tempdir, files{k}]);
 %       end
    
    end

