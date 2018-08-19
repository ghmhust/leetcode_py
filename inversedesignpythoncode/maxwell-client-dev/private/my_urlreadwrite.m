%% maxwell_urlreadwrite
% Initiates an http or https connection to a server.

function [urlConnection, errorid, errormsg] = my_urlreadwrite(urlChar, varargin)
    % Default output arguments.
    urlConnection = [];
    errorid = '';
    errormsg = '';

    % Determine the protocol (before the ":").
    protocol = urlChar(1:min(find(urlChar==':'))-1);

    % Try to use the native handler, not the ice.* classes.
    use_maxwell_cert = false;
    switch protocol
        case 'http'
            try
                handler = sun.net.www.protocol.http.Handler;
            catch exception %#ok
                handler = [];
            end
        case 'https'
            handler = sun.net.www.protocol.https.Handler;
        otherwise
            handler = [];
    end

    % Create the URL object.
    try
        if isempty(handler)
            url = java.net.URL(urlChar);
        else
            url = java.net.URL([],urlChar,handler);
        end
    catch exception %#ok
        errorid = ['InvalidUrl'];
        errormsg = 'Either this URL could not be parsed or the protocol is not supported.';
        return
    end

    % Get the proxy information using MathWorks facilities for unified proxy
    % prefence settings.
    mwtcp = com.mathworks.net.transport.MWTransportClientPropertiesFactory.create();
    proxy = mwtcp.getProxy(); 


    % Open a connection to the URL.
    if isempty(proxy)
        urlConnection = url.openConnection;
    else
        urlConnection = url.openConnection(proxy);
    end

    %% Set timeout.
    % 15-second window to establish connection.
    urlConnection.setConnectTimeout(15e3);    
    % 15-minute window to read from connection.
    urlConnection.setReadTimeout(15e3); 
end

