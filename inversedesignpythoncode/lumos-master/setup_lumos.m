%% setup_lumos
% A script that gets things going for lumos.

    %% Add all subdirectories to path.
    lumos_root_dir = strrep(mfilename('fullpath'), '/setup_lumos', '');
    path(path, genpath(lumos_root_dir));
    % str = strrep(str1, str2, str3)����str1���ҵ�str2���滻��str3
    % p=genpath('a')%��ȡ�ļ���a��·���ַ������Լ��ļ���a�ĸ������ļ��е�·�������ص�p
    % path(path,'a');%��·��a��ӵ�����Ŀ¼�ĵ׶ˣ����a�Ѿ�����������Ŀ¼�У������ƶ����׶�
    %% Install maxwellFDS (alpha).
    % If this fails - try, try again.
    % try-catch-end�ṹ��try ����ִ�е����E   catch  ���E���д���ִ��catch��end֮��Ĵ��� end
    % ����Ϊ��ʱ��break��ѭ���еĹ�����������ǰѭ����continue�Ĺ����ǽ�������ѭ��������һ��ѭ��
%     for i = 1 : 20
%         try
%            eval(urlread('http://m.lightlabs.co/v1'));
%            % ������ַ������jesse���õ��Ʒ�������ַ��maxwellFDs����ָ��pathonд��maxwell fdfd��
%              % eval(urlread('http://m.lightlabs.co/v1'));
%             % S =urlread('URL','method',PARAMS)��urlread�������Զ�ȡ��ҳ����һ������ҳ��ַ���ڶ�����get����post������������Ҫ����ҳ���ݵĲ���
%             break;
%         catch exception % �쳣
%             
%             fprintf(getReport(exception, 'extended'));
%             fprintf('\n');
%             % fprintf(fid,format,variables)����ָ���ĸ�ʽ��������ֵ�������Ļ��ָ���ļ�����ָ���ĸ�ʽ��������ֵ�������Ļ��ָ���ļ�
%             % format����ָ���������ʱ���õĸ�ʽ�� ���� ��\n���У�
%             pause(rand(1)*10);
%             % pause(n)�����÷����ڼ���ִ��ǰ��ִֹ�г���n��
%             % pause(n)�����÷����ڼ���ִ��ǰ��ִֹ�г���n��
%             continue;
%         end
%     end
