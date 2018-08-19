% Script for online activation of maxwell-client.
% This should be uploaded to S3.
% The user can then perform online activation using some form of
%   
%   eval(urlread('http://this-file-right-here-that-you-uploaded-to-the-web'))
%

version = 'dev';

%% Generic install script.
fprintf('Loading maxwell-client (version %s)...', version);

% Some constants.
maxwelldir = [tempdir, filesep, 'maxwell-client'];

% Make the directory.
try 
    warning off;
    rmdir(maxwelldir, 's');
    warning on;
end
mkdir(maxwelldir);

% Get the zip files.
unzip(['https://codeload.github.com/JesseLu/maxwell-client/zip/', version], maxwelldir);
path(genpath(maxwelldir), path);
fprintf('done\n');
