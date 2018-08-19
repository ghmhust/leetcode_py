% Add the current directory (maxwell-client/private) to Matlab's java class path.
function modify_javapath()
    javaaddpath(strrep(mfilename('fullpath'), mfilename(), ''))
