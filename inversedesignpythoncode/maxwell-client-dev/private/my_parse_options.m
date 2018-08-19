function [options] = my_parse_options(options, args, fname)

    for k = 2 : 2 : length(args)
        if isfield(options, args{k-1})
            % isfield�ж������Ƿ��ǽṹ��������򣨳�Ա��
            options = setfield(options, args{k-1}, args{k});
        else
            error('''%s'' is not a valid optional argument for %s', ...
                    args{k-1}, fname);
        end
    end
