function [arrdef,edims,pdims,Idx] = readarraydefinition(filename)
% [ARRDEF,EDIMS,PDIMS,IDX] = READARRAYDEFINITION(FILENAME)
% Wrapper function to import one of three (provisionally four) formats of array definition files:
%
%   A. PREPROCESSOR V2.0 specification: Delimited file with format:
% 
%       #file_format: PREPROCESSOR_V2.0
%       #arrdef_id: some_id
%       #layout_id: some_other_id
%       #comments
%       #...
%       #inv;box2;box1;str;mount;pos
%       'a';'b';'c';'d';1;1
%       ...
% 
%     inv, box1, box2, and str will be read as string-labels; mount, and pos should be unique 
%     one-based index-pairs. For historical reasons, '/' will be treated as empty, in any column.
% 
%   A(bis). PROVISIONAL: since getting the propper column-order and indexing seems to be just too 
%     much to ask from the Pre-Processor, we're currently handling:
% 
%       inv;box1;box2;mount;str;mod
%       'a';'c';'b';0;'d';0
% 
%     Where mount and mod are 0-based indices, and the module position index is guessed from the 
%     order of the rows.
%
%      PROVISIONAL: format ID criterium is exact headerline: 
%                   'inverter;combiner a;combiner b;mount;string;modul'
%
%   B. Delimited file with (with optional header line), at least 2 columns of data:
%       cell;str - used for module definitions.
%       mod;str;inv - used for tracker-inverters (mount = inv, pos = Nstr·(str-1)+ mod)
%       mod;str;inv;pos;mount; - minimum complete PV-plant definition
%       mod;str;box;inv;pos;mount; 
%       mod;str;box;mppt;inv;pos;mount; 
%
%       mod;str;[box1];[box2];[inv];...;@; pos; mount;[row];[field];... - unlimited dimensions
%               note that an '@' character separates the electrical and physical definitions
%
%   C. XLS/XLSX single page with data table 2 to 7 columns (interpreted as in B.)
%   D. XLS/XLSX 2 to 7 page file, each page with an array of integers of the same size,
%      interpreted as the columns of B. and C.
%
% OUTPUT:
%
% ARRDEF - Ne+Np-column array of indices, each row defining the unique electrical coordinates of a
%   module (e.g. inverter, string, module) followed by its unique physical coordinates (e.g. mount,
%   position).
%
% EDIMS - string of electrical dimension 'types' describing the interpretation of the first Ne rows  
%       of arrdef. Allowed characters (types) are:
%             x: absolute index - must come first, only once - currently ignored!
%             c: cell - currently only expected as part of a module definition: [x]cs
%             s: string (DC series) - currently only expected (once), right after c/m
%             m: module
%             b: box (DC parallel) - currently only expected between m and e/i, i.e. conn. box(es)
%             e: MPPT input - if declared, must be followed by i
%             i: inverter
%             t: Trafo - nested trafo/busbar arrangements are allowed
%             a: AC busbar 
%
%       e.g. 'msi' for module-string-inverter
%
% PDIMS - string of physical dimension 'types' describing the interpretation of the last Np rows of . 
%       arrdef. Allowed characters (types) are:
%             x: absolute index - must come first, only once - currently ignored!
%             p: module position - must be first or follow x
%             t: mount / tracker - must follow p
%             b: block - can follow t, nested blocks allowed.
%       
%       e.g. 'pt' for position-tracker
%
% IDX - LABELMAP object containing the univocal correspondence bewteen string & inverter indices
%       and their original (char) labels.
%
% PROVISIONAL: For preprocessor-based formats, connexion-box coordinates are dumped to IDX, and
%       ARRDEF reduced to standard 5-columns: msipt
%
% FUTURE: check/use column header names (cell,mod,str,inv,trafo,series,parallel...)
%         paired-column reduced defintions (mod-str / str-box / box-inv @ mod-trck, trck-block..)

    % browse for file if necessary
	if nargin < 1, filename = '*'; end
    filename = pickfile(filename,'Pick an array-definition file');

    edims = []; % Number of electrical dimensions (empty = unknown)
    pdims = [];
    
    if isxlsarraydef(filename)
    % Call readxlsarraydef for any excel files
        arrdef = readxlsarraydef(filename,false);
    elseif ispparraydef(filename)
    % try with Pre-Processor output [i b1 b2 t ~.~.s m]
        [arrdef,Idx,edims,pdims] = readpparraydef(filename);
    else
    % otherwise assume standard format [m s i p t b1 b2] / [e1,e2...@ p1, p2...]
        [arrdef,edims,pdims] = readsplitarraydef(filename);
        Idx = [];
    end
end

function itis = isxlsarraydef(filename)
% Try to guess and read file type by extension
    fileext = regexp(filename,'\.\w*$','match');
    itis = any(strcmp(fileext{1},{'.xls','.xlsx','.xlsm','.xlsb'}));
end

function itis = ispparraydef(filename)
% PROVISIONAL criteria is that the header reads 'inverter;combiner a;combiner b;mount;string;modul'
% with minor deviations (e.g. spaces, combiner_n, etc.)

	fileID = fopen(filename);
	if fileID < 0, error('readpparraydef:fopen','readpparraydef: Could not open file'); end
	firstline = textscan(fileID,'%[^\n\r]',1);
    fclose(fileID);
    
    expr = '^inverter[;\s]+combiner.*[;\s]+combiner.*[;\s]+mount[;\s]+string[;\s]+modul[e]?$';    
    %itis = strcmp(firstline{1},'inverter;combiner a;combiner b;mount;string;modul');
    itis = ~isempty(regexpi(firstline{1},expr));
end

function [arrdef,edims,pdims] = readsplitarraydef(filename)
% Read a delimited file (with optional header) and format:
%
%   mod;str;[box1];[box2];[inv];...;@; pos; mount;[row];[field];... - unlimited dimensions
%   note that an '@' character separates the electrical and physical definitions
%
% Return an integer array of indices, and the number of electrical and physical dimensions 
% (columns before and after @). If separator '@' is not included, edims = pdims = []
%
% FUTURE: try to infere dimension meanings from headers, and return edims, pdims as strings with
%         code letters (see checkedims / checkpdims).

    % Try to open the file and read the first two lines
    fileID = fopen(filename);
    assert(fileID >= 0, 'readsplitarraydef:fopen','Could not open file');
    foo = onCleanup(@() fclose(fileID));
    firstlines = textscan(fileID,'%[^\n\r]',2);
    delete foo;  % fclose(fileID);

    % Guess delimiters from second line
    delim = [char(9),';, ']; % delimiter priority: (tab) ; , (space)
    for j = 1:numel(delim)
        if any(firstlines{1}{2}==delim(j))
            delim = delim(j);
            break
        end
    end

    % Try to split and read a number
    fields = strsplit(firstlines{1}{1},delim);
    skipheader = isnan(str2double(fields{1}));
    
    % Try to find separator '@'
    fields = strsplit(firstlines{1}{2},delim);
    if firstlines{1}{2}(end) == delim, fields = fields(1:end-1); end  % remove empty field when row closes with delim
    isat = strcmp('@',fields); 
    if ~any(isat)
    % If separator '@' is not found, call dlmread and hope for the best
        arrdef = dlmread(filename,delim,skipheader*1,0);
        edims = []; pdims = [];
    else
        if sum(isat) > 1, error('don''t know what to do with more than one @'); end
        
        formatstring = strjoin(repmat({'%d'},1,numel(fields)),delim);

        fileID = fopen(filename);
        assert(fileID >= 0, 'readsplitarraydef:fopen','Could not open file');
        foo = onCleanup(@() fclose(fileID));
        data = textscan(fileID,formatstring,'Delimiter',delim,'treatasempty',{'/','@'}, ...
                'CollectOutput',1,'HeaderLines',1*skipheader);

        arrdef = data{1}(:,~isat);
        edims = find(isat)-1;
        pdims = size(arrdef,2)-edims;
    end
end

function [arrdef, StrIdx, edims, pdims] = readpparraydef(filename)
% Import an array-definition file as defined in the PREPROCESSOR V2.0 specification:
%
%   #file_format: PREPROCESSOR_V2.0
%   #arrdef_id: some_id
%   #layout_id: some_other_id
%   #comments
%   #...
%   #inv;box2;box1;str;mount;pos
%   'a';'b';'c';'d';1;1
%   ...
%
% inv, box1, box2, and str will be read as string-labels; mount, and pos should be unique one-based
% index-pairs. For historical reasons, '/' will be treated as empty, in any column.
%
% Alternatively, since getting the column-order and indexing right seems to be just too much
% for the people working in the Pre-Processor, we're currently handling:
%
%   inv;box1;box2;mount;str;mod
%   'a';'c';'b';0;'d';0
%
% Where mount and mod are 0-based indices, and the module position index is guessed from the order
% of the rows.
%
% Returns an NDmap object, with the connection scheme of the PV-plant

    % Attempt to read the cannonical format:
    filedata = readtxtfile(filename,{'file_format','arrdef_id','layout_id'},...
        'formatstring','%s %s %s %s %d %d','TreatAsEmpty','/','EmptyValue',NaN,'CollectOutput',0);
    
    trueppformat = isfield(filedata.params,'file_format') && ...
        strcmpi(filedata.params.file_format, 'PREPROCESSOR_V2.0');

    if ~trueppformat
    % PROVISIONAL: fix column order and zero-based indexing
        filedata = readtxtfile(filename,'formatstring','%s %s %s %d %s %d','TreatAsEmpty','/',...
            'EmptyValue',NaN,'CollectOutput',0);
        filedata.data = filedata.data([1 3 2 5 4 6]); % fix column order
        filedata.data{5} = filedata.data{5} + 1;      % fix 0-based mount & module indices 
        filedata.data{6} = filedata.data{6} + 1;
    end
    filedata = filedata.data;
 
    Nm = numel(filedata{1});
    labeleddims = 1:4;
    numbereddims = 5:6;
    
    emptydims = false(Nm,6);
    for j = labeleddims
       emptydims(:,j) = strcmp(filedata{j},'/') | arrayfun(@(x) isempty(x{1}),filedata{j});
       filedata{j}(emptydims(:,j)) = {''};
    end
    
    if ~trueppformat
    % PROVISIONAL: Get module-position number inside mount
        [pcoords,sortidx] = sortrows([filedata{5:6}]);
        [~,ia] = unique(pcoords(:,1),'stable');         % index of first occurence of each mount
        Np = max(int32(diff([ia;Nm+1])));               % modules per mount, based on ia
        if any(pcoords(:,2)-(pcoords(:,1)-1)*Np > Np)
            error('readpparraydef:Npcross','readpparraydef: unexpected module numbering convention.')
        else
            filedata{6}(sortidx) = mod(pcoords(:,2)-1,max(Np))+1;
        end
    end
    
    nanornegative = @(x) isnan(x) | (x < 0);
    emptydims(:,numbereddims) = nanornegative([filedata{numbereddims}]);
    
    % Remove disconnected rows
    notconn = any(emptydims(:,[1,4]) > 0, 2); 
    if any(notconn)
        Nm = Nm - nnz(notconn);
        for j = 1:6, filedata{j} = filedata{j}(~notconn,:); end
        warning('readpparraydef:conn','%d disconnected elements (%s) found and removed',...
                    nnz(notconn),shortliststr(find(notconn)+1,'row',5));
        if all(notconn)
            error('readpparraydef:noconn',...
                '... actually, no modules seem connected, data might have the wrong format');
        end
    end
    
    % Modules that should be connected but have no physical-address are a different thing...
    notthere = any(emptydims(:,5:6) > 0, 2);
    if any(notthere)
       error('readpparraydef:notthere','%d ''connected'' elements (%s) lead nowhere',...
                nnz(notthere),shortliststr(find(notthere)+1,'row',5));
    end
    
    [idx.t,idx.p] = deal(filedata{5:6});
    
    % Get inverter index
    [~,~,idx.i] = unique(filedata{1},'sorted');
    
    % Get minimal string indices
    idx.s = zeros(numel(idx.i),1);
    strlabels = [filedata{2:4}];
    for i = 1:max(idx.i)
        [~,~,idx.s(idx.i == i)] = uniquecell(strlabels(idx.i == i,:),'rows','sorted');
    end
    
    % Create a combined inverter-string label index
    [~,ia,ic] = unique([idx.i,idx.s],'rows','sorted');
    
    % ... Remove empty dimensions first
    labeleddims = nonzeros(~all(emptydims(:,labeleddims),1).*labeleddims);
    strlabels = [filedata{labeleddims}];
    StrIdx = LabelMap([idx.i(ia),idx.s(ia)],strlabels(ia,:));
    
    % Quick check number of modules per string in each inverter
    nm = accumarray(ic,1);
    errmsg = '';
    badmodules = false(numel(ic),1);
    for j = 1:StrIdx.esize(1)
        f = idx.i(ia) == j;
        bad = nm(f) ~= mode(nm(f));
        if any(bad)
            badmodules = badmodules | idx.i == j & any(idx.s == find(bad)',2);
            c = StrIdx.plabels(j,bad);
            errmsg = sprintf('%s\n\tInverter %s (%s modules): Check %s',...
                errmsg,c{1},shortliststr(nm(f),'and','/'),shortliststr(c(:,2),'string'));
        end
    end
    if ~isempty(errmsg)
        try %#ok TRYNC
            badmodules = double(unique([idx.t(badmodules),idx.p(badmodules)],'rows'));
            Trck = evalin('base','Trackers');
            bad = full(sparse(badmodules(:,1),badmodules(:,2),1,size(Trck.centers,2),Trck.geom.dims(1)));
            GUIfigure('error','Possible Errors');
            plotarrayprop(Trck,bad);
        end
        error('Inconsistent number of modules withing strings of the same inverter:%s',errmsg);
    end
    
    % Get module number inside string
    arrdef = sortrows([idx.i,idx.s,idx.t,idx.p]); % Sort by inverter-string-mount-position
    [~,idx,ida] = unique(arrdef(:,1:2),'rows','stable');
    arrdef(:,5) = int32((1:Nm)' - idx(ida) + 1);
    
    % Arrange according to (deprecated) ArrayDefinition constructor convention: [m,s,i,p,t] 
    arrdef = arrdef(:,[5 2 1 4 3]);
    
    % if everything worked, return dimension code letters
    edims = 'msi';
    pdims = 'pt';
end

function arraydef = readxlsarraydef(filename,originalsize)
% Reads a k-depth array definition (i.e. a set of 1 to k integer indices for every
% element of an array, which can be used to group these in a hierarchical structure). 
% indices are read from an excel file with one to k sheets, each ith an m·n matrix of integers.
%
% Elements for which the first sheet contains zero will be ignored on all other sheets.
%
% Returns an (mn)·k array where each row should be a unique vector [element,group,super-group,...],
% and rows are sorted according to the main index.
% i.e. arraydef(j,h) represents the index of element find(index==j) at hierarchy level h
%
% If originalsize == true, the output is reshaped into an (m,n,k) array, that is easier to verify

    if nargin < 2, originalsize = false; end

	% Read file, and get the number of non-empty pages 
	[status,sheets] = xlsfinfo(filename);
    if isempty(status)
		error('Cannot read specified file');
    end
    k = numel(sheets);
	while isempty(xlsread(filename,sheets{k},'','basic'))
		k = k-1; 
	end
	if k < 1, error('File appears empty, or has no numeric data'); end
	
	% Read main-index, and check it
	data = xlsread(filename,1,'','basic');		% Read first page
    
    % Read single-page files as a table
    if k == 1
        if size(data,2) < 3 || size(data,2) > 8
            error('Definition table doesn''t have a recognizable format')
        else
            arraydef = data;
            return;
        end
    end
    
	[N,M] = size(data);
	arraydef = zeros(N*M,k);		
	arraydef(:,1) = data(:);

	% Read any additional planes
    for j = 2:k
		data = xlsread(filename,sheets{j},'','basic');
		if any(size(data) ~= [N,M])
			error('readarraydefinition:arraysize','Array size mismatch')
		else
			arraydef(:,j) = data(:);
		end
    end
    
    used = all((arraydef > 0)& ~isnan(arraydef),2);
	arraydef = arraydef(used,:);
    
	% DEBUG: for easy comparison
	if originalsize
		newarraydef = zeros(N*M,k);
		newarraydef(used,:) = arraydef;
		arraydef = reshape(newarraydef,[N,M,k]);
	end
end

		
