function writeoptionsfile(path)
% Write a modified copy of the DefaultOptions.m function as ./AddSimOptions.m on the current
% directory, so that the user can modify it.

    here = pwd(); lastwill = onCleanup(@() cd(here));
    if nargin > 0 && ~isempty(path), cd(path); end
    
    txt = regexp(fileread('DefaultOptions.m'),'\n|\r|\n\r','split')';

    % remove function header/end, and .version, .prjname fields on start
    for a = 1:numel(txt), if ~isempty(strfind(txt{a},'<USER_OPTIONS>')), break; end; end
    for b = numel(txt):-1:1, if ~isempty(strfind(txt{b},'</USER_OPTIONS>')), break; end; end
    txt = txt(a+1:b-1);
    
    txt = strrep(txt,'Defaults.','SimOptions.');

    % Backup exsisting AddSimOptions.m file
    if ~isempty(dir('AddSimOptions.m')), copyfile('AddSimOptions.m','AddSimOptions.m~'); end

    fID = fopen('AddSimOptions.m','w');
    lastwill(2) = onCleanup(@() fclose(fID)); %#ok<NASGU>
    
    fprintf(fID,'\n\tglobal SimOptions;\n');
    cellfun(@(x) fprintf(fID,'%% %s\n',x),txt);
end
