function filename = writeoptionsfile(filename,overwrite)
% WRITEOPTIONSFILE(PATH) - print global SimOptions structure to PATH/simoptions.json

    narginchk(0,2);
    if nargin < 2, overwrite = false; end

    global SimOptions;
    if isempty(SimOptions), SimOptions = completeoptions(); end

    if nargin < 1 || isempty(filename)
        [path,name] = fileparts(SimOptions.prjname);
        if isempty(name)
            name = 'simoptions.json'; 
        else
            name = [name '_simoptions.json']; 
        end
        filename = fullfile(path,name);
    end
    
    try
        txt = jsonencode(SimOptions);
    catch ERR
        if ~contains(ERR.message,'function_handle'), rethrow(ERR); end
        [val,fld] = nestedstruct2cell(SimOptions);
        isfcn = cellfun(@(x) isa(x,'function_handle'),val);
        if ~any(isfcn), rethrow(ERR); end
        for j = find(isfcn)'
            SimOptions = setnestedfield(SimOptions,fld{j},func2str(val{j}));
        end
        txt = jsonencode(SimOptions);
    end
    txt = regexprep(txt,'("\w*":.*?)(?="\w*":)',['$1',newline()]);

    if isfile(filename) && ~overwrite, backupdelete(filename,'-warning'); end

    fID = fopen(filename,'w');
    lastwill = onCleanup(@() fclose(fID));

    fprintf(fID,'%s\n',txt);
end
