%% setup_lumos
% A script that gets things going for lumos.

    %% Add all subdirectories to path.
    lumos_root_dir = strrep(mfilename('fullpath'), '/setup_lumos', '');
    path(path, genpath(lumos_root_dir));
    % str = strrep(str1, str2, str3)，在str1中找到str2，替换成str3
    % p=genpath('a')%获取文件夹a的路径字符串，以及文件夹a的各级子文件夹的路径，返回到p
    % path(path,'a');%将路径a添加到搜索目录的底端，如果a已经存在于搜索目录中，则将其移动到底端
    %% Install maxwellFDS (alpha).
    % If this fails - try, try again.
    % try-catch-end结构：try 尝试执行的语句E   catch  如果E运行错误，执行catch和end之间的代码 end
    % 条件为真时，break在循环中的功能是跳出当前循环，continue的功能是结束本次循环跳到下一次循环
%     for i = 1 : 20
%         try
%            eval(urlread('http://m.lightlabs.co/v1'));
%            % 这里网址可能是jesse所用的云服务器网址，maxwellFDs可能指代pathon写的maxwell fdfd包
%              % eval(urlread('http://m.lightlabs.co/v1'));
%             % S =urlread('URL','method',PARAMS)，urlread函数可以读取网页，第一个是网页地址，第二个是get或是post，第三个则是要向网页传递的参数
%             break;
%         catch exception % 异常
%             
%             fprintf(getReport(exception, 'extended'));
%             fprintf('\n');
%             % fprintf(fid,format,variables)，按指定的格式将变量的值输出到屏幕或指定文件，按指定的格式将变量的值输出到屏幕或指定文件
%             % format用来指定数据输出时采用的格式， 换行 （\n换行）
%             pause(rand(1)*10);
%             % pause(n)：此用法将在继续执行前中止执行程序n秒
%             % pause(n)：此用法将在继续执行前中止执行程序n秒
%             continue;
%         end
%     end
