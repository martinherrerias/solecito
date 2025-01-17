function varargout = readoptionsfile(file)
% OPT = READOPTIONSFILE(FILE) - read options FILE (as generated by WRITEOPTIONSFILE.
% READOPTIONSFILE(FILE) - read options FILE and apply to global SimOptions

    if ~isfile(file), error('Failed to find file: %s',file); end
    txt = fileread(file);
    txt = regexprep(txt,'#[^\n\r]*',''); % remove comments
    opt = jsondecode(txt);

    [val,fld] = nestedstruct2cell(opt);
    isfcn = cellfun(@(x) (ischar(x) || isstring(x)) && numel(x) > 2 && strcmp(x(1:2),'@('),val);
    if any(isfcn)
        for j = find(isfcn)'
            opt = setnestedfield(opt,fld{j},str2func(val{j}));
        end
    end
    isvec = cellfun(@(x) isvector(x) && size(x,1) > 1,val);
    if any(isvec)
        for j = find(isvec)'
            opt = setnestedfield(opt,fld{j},val{j}');
        end
    end

    if nargout == 0
        completeoptions(opt);
    else
        varargout = {opt};
    end
end
