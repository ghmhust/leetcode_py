% Batch download via http.
function my_download(filenames, local_dir, url)
    for k = 1 : length(filenames)
        [inputStream{k}, file{k}, conn{k}] = my_open_connection(filenames{k}, local_dir, url);
    end

    my_stream_send(inputStream, file);

    for k = 1 : length(filenames) % Close the files and connections.
        inputStream{k}.close();
        file{k}.close();
        conn{k}.disconnect();
    end

end


function [inputStream, file, urlConnection] = my_open_connection(filename, local_dir, url)
    url = [url, filename];
    urlConnection = my_urlreadwrite(url);
    urlConnection.connect();
    inputStream = urlConnection.getInputStream();
    file = java.io.FileOutputStream([local_dir, filesep, filename]);
end
                 

