function my_update_documentation()
% Update documentation for Maxwell.
    
    maxwell_dir = strrep(mfilename('fullpath'), ['private/', mfilename], '');
    files = dir(maxwell_dir);

        %
        % Find example_*.m and maxwell_*.m files
        %

    [maxwell_files, example_files, maxopt_case_files, maxopt_example_files, maxopt_files] = deal({});
    for k = 1 : length(files)
        if strfind(files(k).name, 'maxwell_') == 1 & ...
                files(k).name(end-1:end) == '.m' & ...
                ~strcmp(files(k).name, 'maxwell_help.m')
            maxwell_files{end+1} = files(k).name;

        elseif strfind(files(k).name, 'example') == 1 & ...
                files(k).name(end-1:end) == '.m'
            example_files{end+1} = files(k).name;

        elseif strfind(files(k).name, 'maxopt_case') == 1 & ...
                files(k).name(end-1:end) == '.m'
            maxopt_case_files{end+1} = files(k).name;

        elseif strfind(files(k).name, 'maxopt_example') == 1 & ...
                files(k).name(end-1:end) == '.m'
            maxopt_example_files{end+1} = files(k).name;

        elseif strfind(files(k).name, 'maxopt') == 1 & ...
                files(k).name(end-1:end) == '.m'
            maxopt_files{end+1} = files(k).name;
        end
    end
    doc_files = {maxwell_files, example_files, ...
                maxopt_files, maxopt_example_files, maxopt_case_files};


        %
        % Generate maxwell_help file.
        %

    % Get template.
    f = fopen('my_help_template.m', 'r');
    my_help = char(fread(f))';
    fclose(f);

    section_name = {'<<<functions>>>', '<<<examples>>>', ...
                    '<<<maxopt_functions>>>', '<<<maxopt_examples>>>', ...
                    '<<<maxopt_cases>>>'};

    % Get by-line (second line) from files.
    for j = 1 : length(doc_files)
        oneliners = {};
        for k = 1 : length(doc_files{j})
            f = fopen([maxwell_dir, doc_files{j}{k}], 'r');
            title = fgetl(f); % Usually not needed.
            byline{k} = fgetl(f);
            try
                byline{k}(1:2) = []; % Delete first 2 characters (should be '% ').
            end 
            fclose(f);
            % Make one-liner.
            oneliners{k} = sprintf('%%\n%% * <%s %s> - %s\n', ...
                                    strrep(doc_files{j}{k}, '.m', '.html'), ...
                                    strrep(doc_files{j}{k}, '.m', ''), ...
                                    byline{k});
        end
        my_help = strrep(my_help, section_name{j}, [oneliners{:}]);
    end

    f = fopen([maxwell_dir, 'maxwell_help.m'], 'w');
    fwrite(f, my_help);
    fclose(f);


        %
        % Publish documentation.
        %

    for j = 1 : length(doc_files)
        for k = 1 : length(doc_files{j})
            publish([maxwell_dir, doc_files{j}{k}], 'evalCode', false);
        end
    end

    publish([maxwell_dir, 'maxwell_help.m'], 'showCode', false, 'evalCode', false);

