function [status,msginfo] = writearraydefinition(obj,filename,header,forceoverwrite)
% [STATUS,MSGINFO] = WRITEARRAYDEFINITION(OBJ,FILENAME,HEADER,FORCEOVERWRITE)
% Writes an array definition file (csv table or multi-page xls file) describing the structure
% of a pvArea-object, NDmap object, or integer array obj. 
% The format of the output is given by the filename extension:
%
%   .xls* - Multi-page file
%   .csv - Multiple files, one for each plane / depth level
%   .arrdef - Pre-processor V2.0 convention: #header + [i;b2;b1;s;t;p]
%   others - csv table. If no extension is provided, .arrdef is used.
% 
% If OBJ is a pvArea or NDmap object, the output will be such that 
%   OBJ = reshapearea(flattenpvarea(OBJ),readarraydefinition(filename)).
%   Each column/sheet of the produced file will contain a set of integers, representing, in order:
% 	1. index - unique numbers for each element (sort order, for reference to element list of flattenpvarea(obj))
% 	2. element - unique indices inside each group
% 	3. group - repeated integers for all elements which belong to the same group.
% 	4. super-group - repeated integers for all elements which belong to the same super-group...
% 	...
% 	k+1. largest hierarchical group
%
% If OBJ is an N·M·k | k > 1, integer array, sheets/columns will be created for each plane (:,:,j)
% If OBJ is an NM·k integer array, it will be printed as such in a single table.
%
% For .arrdef output, input rows/planes are expected in standard order [m s i p t b1 b2]
% with at least the first 5 defined (b1 and b2 set to 1 if missing). ';' is used as delimiter.
% HEADER in this case should be left empy, or be a structure with fields arrdef_id and layout_id.
%
% For all other formats, a string HEADER should be provided, to set the delimiter character,
% propperly label columns/sheets, and separate electrical and physical coordinates.
% Ideally, the HEADER should have the form:
%
%       'mod;str;inv;box... ;@;pos;trck;blck...'
%
%       a. The number delimited fields shoud match the number of columns/planes in the array
%       b. The optional '@' character (not counted as field) separates electrical and physical
%          dimensions.
%       c. Default header for pvAreas is 'el;gr;sgr;ssgr;sssgr;...'
%          Default header for NDmap objects is 
%              'mod;str;[box];[sbox]...;[mppt];inv;@;pos;trck;[blck];[sblck]...'
%       c. Default headers for integer arrays (depending on number of columns/planes) are:
%           2. 'cell;str'                
%           5. 'mod;str;inv;@;pos;mount;
%           7. 'mod;str;box;box;inv;@;pos;mount;
%
% FUTURE: n-dimensional physical & electrical array definitions (e.g. e1,e2,e3,e4..@,p1,p2,p3..)
%         paired-column reduced defintions (mod-str / str-box / box-inv @ mod-trck, trck-block..)
    
    % Set default delimiter priority: (tab) ; , (space)
    try DELIMLIST = getSimOption('DelimiterPriority'); catch, DELIMLIST = [char(9),';, ']; end

    if nargin < 2 || isempty(filename), filename = datestr(now,'yymmdd_HHMM_SSFFF'); end
    if nargin < 3, header = ''; end 
    if nargin < 4, forceoverwrite = false; end
    
    switch class(obj)
    case 'pvArea'
        % Get the array definition, try to fit into pretty-planes
        arrdef = getarraydefinition(obj,true);
    case 'NDmap'
        % Get the generating-table for an already existing NDmap object
        arrdef = deftable(obj);
    otherwise
        if isnumeric(obj)
            arrdef = obj;
        else
            error('Expected pvArea, NDmap, or numeric array');
        end
    end
    
    % Get file type by extension: 'X' = excel,'M' = multi CSV,'P' = preprocesor,'T' = text
    [filename,filetype] = getfiletype(filename,size(arrdef,3));
        
    % Set k as the number of planes/columns
    if any(filetype == ['X','M'])
        k = size(arrdef,3);
    else
        arrdef = flattenarray(arrdef);
        k = size(arrdef,2);
    end
    
    if filetype~='P'
        if isempty(header), header = defaultheader(obj); end
        if iscell(header), header = strjoin(header,DELIMLIST(1)); end
        [colheads,issplit,delim] = workoutheader(header,k,DELIMLIST);
    end
    
    if any(filetype == ['X','M'])
        if isempty(colheads)
            sheetnames = num2str((1:k)','%E02d');
        end
        if filetype == 'M'
            sheetnames = [repmat(filename(1:end-4),k,1),sheetnames,repmat('.csv',k,1)];
        end      
        sheetnames = cellstr(sheetnames);
    elseif filetype ~= 'P' && any(issplit)
            % Put a gap with int32('@') in tabular arrdef
            left = arrdef(:,1:find(issplit)-1);
            right = arrdef(:,find(issplit):end);
            ats = repmat(int32('@'),size(arrdef,1),1);
            arrdef = [left,ats,right];
    end
    
    if filetype ~= 'M' 
        [ok2write,~,msginfo] = right2overwrite(filename,forceoverwrite);
    else
        [ok2write,msginfo] = right2overwrite(sheetnames,forceoverwrite);
    end
    if ~ok2write, status = 0; return; end
            
    switch filetype
        case {'X','x'}
            [status,msginfo] = xlswrite(filename,{''}); % Create file if it doesn't exist
        if filetype == 'x'
            % Write single-plane table in one sheet, update status and message within writemore()
            cleanxlsfile(filename,'arrdef');    
            xlswritemore(filename,colheads,'arrdef');
            xlswritemore(filename,arrdef,'arrdef',xlrange(arrdef,[size(colheads,1),0]));
            xlswritemore(filename,ats,'arrdef',xlrange(ats,[size(colheads,1),find(issplit)]));
        else
            % Write each plane in a new sheet, update status and message within writemore()
            cleanxlsfile(filename,sheetnames);
            sheetnames = sheetnames(~issplit); % Leave @ sheet empty, if existing
            for j = 1:k
                xlswritemore(filename,arrdef(:,:,j),sheetnames{j});
            end
        end
    case 'M'
        % if k ~= 1, filetype would be 'T'
        for j = 1:k
            % Check existance of each page-file
            sheetnames{j} = sprintf('%s%s.csv',filename(1:end-4),sheetnames{j});
            
            if ~isempty(info)
                dlmwrite(sheetnames{j},info,'delimiter',delim,'newline','pc');
                dlmwrite(sheetnames{j},arrdef(:,:,j),'-append','delimiter',delim,'newline','pc');
            else
                dlmwrite(sheetnames{j},arrdef(:,:,j),'delimiter',delim,'newline','pc');
            end
        end
    case 'T'
        
        if ~any(issplit)
            formatstring = strjoin(repmat({'%d'},1,k),delim);
        else
            formatstring = strjoin(repmat({'%d'},1,size(left,2)),delim);
            formatstring = [formatstring,delim,'%c',delim]; %switch field to character
            formatstring = [formatstring,strjoin(repmat({'%d'},1,size(right,2)),delim)];
        end
        formatstring = [formatstring '\r\n'];
        
        fID = fopen(filename,'w');
        if ~isempty(header), fprintf(fID,'%s\r\n',header);end
        fprintf(fID,formatstring,arrdef');
        fclose(fID);
        
    case 'P'
        % Complete array [m s i p t b1 b2]
        [NM,k] = size(arrdef);
        if k < 5, error('Cannot generate pre-processor file for incomplete array definitions'); end
        if k < 6, arrdef = [arrdef ones(NM,2)]; end
        if k < 7, arrdef = [arrdef ones(NM,1)]; end
        
        % Change column order: [m s i p t b1 b2] -> [i b2 b1 s t p]
        arrdef = arrdef(:,[3 7 6 2 5 4]);
        
        %PROVISIONAL! preprocessor not up to specs
        assert(isempty(header) || isstruct(header),'Structure HEADER required for *.arrdef output');
        if isempty(header), header = struct(); end
        arrdef_id = regexprep(filename,'\.[^.]{0,10}$','');      % remove extension
        if isfield(header,'layout_id'), layout_id = header.layout_id; else, layout_id = arrdef_id; end
        
        if isfield(header,'comments')
            comments = cellfun(@(r) sprintf('#  %s',r),header.comments,'unif',0);
        else
            comments = {};
        end
        header = {'#file_format: PREPROCESSOR_V2.0'};
        header{2} = ['#arrdef_id: ' arrdef_id];
        header{3} = ['#layout_id: ' layout_id];
        header = cat(1,header(:),comments(:));

        %header{4} = '#inverter_model_idx: model_1{A;B;C}; model_2{D;E}; ...';
        header{end+1} = '#inverter;combiner a;combiner b;string;mount;position';
        
        fID = fopen(filename,'w');
        cellfun(@(c) fprintf(fID,'%s\r\n',c),header);
        fprintf(fID,'%d;%d;%d;%d;%d;%d\r\n',arrdef'); 
        fclose(fID);
    end

    function success = xlswritemore(file,data,varargin)
        [success,msginfo] = xlswrite(file,data,varargin{:});
        status = status && success;
    end
end

function arrdef = flattenarray(arrdef)
% flatten array definition, if required
    if size(arrdef,3) == 1            
        return;
    else
        [N,M,k] = size(arrdef);
        arrdef = shiftdim(arrdef,2);          % shift to a k·N·M array
        arrdef = reshape(arrdef,[k,N*M])';      % turn into a table, N·M rows by k colums
        arrdef = sortrows(arrdef,k:-1:1);
    end
end

function [filename,type] = getfiletype(filename,depth)
    % Get file extension, if specified, use .arrdef by default
    fileext = regexp(filename,'\.\w*$','match');
    if isempty(fileext)
        filename = [filename '.arrdef']; 
        fileext{1} = '.arrdef';
    end

    switch fileext{1}
    case '.arrdef'
        type = 'P';
    case {'.xls','.xlsx','.xlsm','.xlsb'}
        if depth > 1, type = 'X';
        else, type = 'x';
        end
    case '.csv'
        if depth > 1, type = 'M'; % M for Multiple sheets
        else, type = 'T';                 % T for single text sheet
        end
    otherwise
        type = 'T';                      % T for single text sheet
    end
end
   
function [colheaders,issplit,D] = workoutheader(header,k,DELIMLIST)
% Try to find delimiter and individual column/page headers from string:
% Return a format-string that works with fprintf(formatstr,1:k)
% and a k-cell of sheetnames for use in multi-page files

    D = DELIMLIST(1);
    colheaders = {};
    issplit = false;

    % Find delimiter that makes sense
    foundit = false;
    for j = 1:numel(DELIMLIST)
        D = DELIMLIST(j);
        if any(header == D)
            colheaders = strsplit(header,D);
            colheaders = colheaders(~strcmp('',colheaders))'; % Remove empty
            
            % If number of fields matches k, found it!
            issplit = strcmp('@',colheaders);
            if numel(colheaders) - 1*any(issplit) == k
                foundit = true;
                break;
            end
        end
    end
    
    if ~foundit % Dump header with a warning
        warning('writearraydefinition:header','header doesn''t match array');
    end
end

function header = defaultheader(obj)

    switch class(obj)
    case 'pvArea'
    % 'el;gr;sgr;ssgr;sssgr;...'
        header{1} = 'el';
        for j = (obj.depth-1):-1:2
            header{j} = sprintf('%sgr',repmat('s',1,j-2));
        end
    case 'NDmap'
    % 'mod;str;[box];[sbox]...;[mppt];inv;@;pos;trck;[blck];[sblck]...'
        e = numel(obj.esize);
        if e > 0, ehead{1} = 'mod'; end
        if e > 1, ehead{2} = 'str'; end
        if e > 2, ehead{e} = 'inv'; end
        if e > 3, ehead{e-1} = 'mppt'; end
        for j = (e-2):-1:3
            ehead{j} = sprintf('%sbox',repmat('s',1,j-3));
        end

        p = numel(obj.psize);
        if p > 0, phead{1} = 'pos'; end
        if p > 1, phead{2} = 'trck'; end
        for j = p:-1:2
            phead{j} = sprintf('%sblck',repmat('s',1,j-2));
        end
        
        header = {ehead{:},{'@'};phead{:}};
        
        otherwise
        % Only special cases
        	k = size(obj,3);
            if k == 1, k = size(obj,2); end
            
            switch k
                case 2, header = {'cell','str'};
                case 5, header = {'mod','str','inv','@','pos','mount'};
                case 7, header = {'mod','str','box','box','inv','@','pos','mount'};
                otherwise
                    header = {};
            end
    end
end
    
function txt = xlrange(obj,offset)
% Return the excel range address (A1 format) of an array obj, starting at position offset (r,c).

    if nargin < 2, offset = [0,0]; end
    
    r = [1,size(obj,1)] + offset(1);
    c = [1,size(obj,2)] + offset(2);
    
    txt = sprintf('%s%d:%s%d',base27dec(c(1)),r(1),base27dec(c(2)),r(2));
end

function d = base27dec(s)
% Taken from XLSWRITE
%   BASE27DEC(S) returns the decimal of string S which represents a number in
%   base 27, expressed as 'A'..'Z', 'AA','AB'...'AZ', and so on. Note, there is
%   no zero so strictly we have hybrid base26, base27 number system.
%
%   Examples
%       base27dec('A') returns 1
%       base27dec('Z') returns 26
%       base27dec('IV') returns 256
%-----------------------------------------------------------------------------

    if length(s) == 1
       d = s(1) -'A' + 1;
    else
        cumulative = 0;
        for i = 1:numel(s)-1
            cumulative = cumulative + 26.^i;
        end
        indices_fliped = 1 + s - 'A';
        indices = fliplr(indices_fliped);
        indices_in_cells = mat2cell(indices, 1, ones(1,numel(indices))); %#ok<MMTC>
        d = cumulative + sub2ind(repmat(26, 1,numel(s)), indices_in_cells{:});
    end
end

function cleanxlsfile(filename,sheets)
% Delete all sheets of Excel file 'filename', leave as an empty file with sheets 'sheets'

    [~, oldsheets] = xlsfinfo(filename);

    objExcel = actxserver('Excel.Application');
    objExcel.DisplayAlerts = false;
    WB = objExcel.Workbooks.Open(fullfile(pwd,filename)); % Full path is necessary!
    
    % Create missing worksheets
    newsheets = setdiff(sheets,oldsheets);
    for i = 1:numel(newsheets)
        lastsheet = WB.WorkSheets.Item(WB.WorkSheets.Count);
        lastsheet = WB.WorkSheets.Add([],lastsheet);
        set(lastsheet,'Name',newsheets{i});
    end

    % Delete existing worksheets that are not in 'sheets'
    sheets2remove = setdiff(oldsheets,sheets);
    for i = 1:numel(sheets2remove)
        WB.Worksheets.Item(sheets2remove{i}).Delete;
    end
    
    % Clear contents of all sheets
    for i = 1:numel(sheets)
        WB.Worksheets.Item(sheets{i}).Cells.Clear;
    end
    
    WB.Save;
    WB.Close;
    objExcel.Quit;
    objExcel.delete;
end


		
