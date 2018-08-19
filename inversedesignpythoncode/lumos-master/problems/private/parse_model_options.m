function [model_options] = parse_model_options(custom_model_options)

    % Default parameters.
    model_options = struct( 'size', 'small', ...
                            'flatten', true, ...
                            'S_type', 'average');

    %% Parse optional parameters.
    for k = 2 : 2 : length(custom_model_options)
        model_options = setfield(model_options, custom_model_options{k-1}, ...
                                                custom_model_options{k});
       % setfield: 给结构数组的字段指定值,s = setfield(s,'field',value)等效于s = setfield(s,'field',value)
    end
end
