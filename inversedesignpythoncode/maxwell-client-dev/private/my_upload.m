% Batch upload via http.
function my_upload(filenames, local_dir, url)

    for k = 1 : length(filenames)
        % Open connections and send headers.
        [infile{k}, outputStream{k}, footer{k}, urlConn{k}] = ...
            my_open_connection(filenames{k}, local_dir, url);
    end

    % Stream the files.
    my_stream_send(infile, outputStream);

    for k = 1 : length(filenames)
        % Send the footers.
        outputStream{k}.write(footer{k}.getBytes(), 0, footer{k}.length); 

        % Close off the messages.
        infile{k}.close(); 
        outputStream{k}.flush();
    end

    % Check for a error response codes.
    for k = 1 : length(filenames)
        if (urlConn{k}.getResponseCode() ~= 200)
            error('%s', char(urlConn{k}.getHeaderField(0)));
        end
    end

    % Close connection.
    for k = 1 : length(filenames)
        outputStream{k}.close();
        urlConn{k}.disconnect();
    end
end

%% Open a connection (POST).
function [infile, outputStream, footer, urlConnection] = ...
            my_open_connection(filename, local_dir, url)
    

    % Create a urlConnection.
    [urlConnection, errorid, errormsg] = my_urlreadwrite(url);
    if isempty(urlConnection)
        error(['Could not connect to url: ', url]);
    end

    urlConnection.setDoOutput(true); % Sets the request mode to POST.
    boundary = '*** maxwell_client boundary ***';
    urlConnection.setRequestProperty('Content-Type',...
        ['multipart/form-data; boundary=', boundary]);

    eol = [char(13),char(10)]; % End-of-line character.

    % Build the header, body and footer of the POST request.
    header = [];
    for k = 1 : 2 : length(params) % Form data for text parameters.
        header = [header, '--', boundary, eol, ...
                    'Content-Disposition: form-data; name="', params{k}, '"', ...
                    eol, eol, params{k+1}, eol];
    end
    header = java.lang.String([header, '--', boundary, eol, ...
                'Content-Disposition: form-data; name="file"; ', ...
                    'filename="', filename, '"', eol, ...
                'Content-Type: application/octet-stream', eol, eol]); 
                % Form data for binary data (the simulation file.
    file = java.io.File([local_dir, filesep, filename]);
    footer = java.lang.String([eol, '--', boundary, '--', eol]);

    % We used a streaming connection, crucial for large files.
    total_length = header.length + file.length() + footer.length;
    urlConnection.setFixedLengthStreamingMode(total_length);
    outputStream = java.io.DataOutputStream(urlConnection.getOutputStream);

    outputStream.write(header.getBytes(), 0, header.length); % Send the header.

    infile = java.io.FileInputStream(file);
end
