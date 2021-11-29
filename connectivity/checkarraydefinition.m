function [arrdef,edims,pdims] = checkarraydefinition(arrdef,edims,pdims,psize)
% [arrdef,edims,pdims] = CHECKARRAYDEFINITION(arrdef)
% [arrdef,edims,pdims] = CHECKARRAYDEFINITION(arrdef,edims,pdims,psize)
% [arrdef,edims,pdims] = CHECKARRAYDEFINITION(arrdef,edims,pdims,pidx)
%
% Take an (Ne+Np)-column arrdef and test to see if it makes a meaningfull array definition, taking
% into account uniqueness, completeness, etc.
%
% edims - string of electrical dimension 'types' describing the interpretation of the first Ne rows  
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
%       If isnumeric(edims), it will be interpreted as the number of electrical dimensions (Ne), and
%       their meaning will be guessed from 'standard' definitions {x, cs, msi, msbi, msbei}
%
% pdims - string of physical dimension 'types' describing the interpretation of the last Np rows of . 
%       arrdef. Allowed characters (types) are:
%             x: absolute index - must come first, only once - currently ignored!
%             p: module position - must be first or follow x
%             t: mount / tracker - must follow p
%             b: block - can follow t, nested blocks allowed.
%       The only recognized 'standard' definitions are currently {x,pt}
%
% NOTE: if neither pdims or edims are known, 'standard' guesses for the Ne+Np columns of arrdef
%       are {cs,msi,msipt,msbipt,msbeipt}
%
% psize - physical layout size, if provided, it will be compared against the Np column indices.
%       IMPORTANT: note that array indices are expected in order child -> parent (e.g. msit, ptb)
%       whereas psize is expected as parent -> child, in agreement with NDmap
%
% pidx - instead of psize (maximum index for each dimension) an extensive list of allowed indices
%       for each physical dimension can be provided as a cell-array of vectors. i.e. pidx{j}
%       contains a list of all valid indices for physical dimension j (column Ne + 1 in arrdef). 
%
% FUTURE: more flexible/complete array definitions (e.g. xcbsmsbbeitat)
%       paired-column reduced defintions (mod-str / str-box / box-inv @ mod-trck, trck-block..)

    if nargin < 4, psize = []; end
    
    % Sort-out incomplete edims / pdims, guess one from the other, or guess
    if nargin < 3 || (nargin > 2 && isempty(pdims)) % unknown pdims ...
        if nargin > 1 && ~isempty(edims) % and known edims
            pdims = size(arrdef,2)-dimsize(edims); 
            assert(pdims >= 0,'checkarrdef:negpdims','Inconsistent arrdef and edims');
        else % both unknown, try to guess
            switch size(arrdef,2)
                case 2, edims = 2; % cell;str
                case 3, edims = 3; % mod;str;inv
                case 5, edims = 3; % mods;str;inv;pos;trck
                case 6, edims = 4; % mod;str;box;inv;pos;trck
                case 7, edims = 5; % mods;str;box;mppt;inv;pos;trck
                otherwise, error('checkarrdef:ncols','Unrecognizable array definition');          
            end
            warning('checkarrdef:nodims','Assuming that %d colums of the array-definition are electrical coordinates',edims);
            pdims = size(arrdef,2)-edims;
        end
    else % known pdims
        if isempty(edims) % ... and unknown edims
            edims = size(arrdef,2)-dimsize(pdims); 
            assert(edims >= 0,'checkarrdef:negedims','Inconsistent arrdef and pdims');
        else % both known, check
            assert(dimsize(edims)+dimsize(pdims) == size(arrdef,2),'Inconsistent arrdef and e+pdims'); 
        end
    end
    
    Ne = dimsize(edims);
    Np = dimsize(pdims);

    % Guess/ask electrical coorinate interpretation from Ne
    [edims, arr] = checkedims(edims,arrdef(:,1:Ne));
    if size(arr,2)~= Ne
        arrdef = [arr arrdef(:,Ne+1:end)]; 
        Ne = size(edims);
    end

    % Guess/ask physical coorinate interpretation from Np
    [pdims, arr] = checkpdims(pdims,arrdef(:,Ne+1:end));
    if size(arr,2)~= Np
        arrdef = [arrdef(:,1:Ne) arr]; 
        %Np = size(arrdef,2)-Ne;
    end
    
    % Remove disconnected
    if Ne > 0
        conn = all(arrdef(:,1:Ne)> 0,2); 
        if any(~conn)
            arrdef = arrdef(conn,:);
            warning('checkarrdef:conn','%d disconnected elements found and removed', sum(~conn));
        end
    end
    
   	% Remove duplicates
    Nm = size(arrdef,1);
    [~,idx] = unique(arrdef,'rows','stable');
    if numel(idx)~= Nm
        warning('checkarrdef:duplicates','%d repeated rows found and removed', numel(idx)-Nm);
        arrdef = arrdef(idx,:);
    end
        
    % Check uniqueness, regularity, dimension order, etc.
    checkelectrical(arrdef(:,1:Ne),edims);
    
    % Check uniqueness, regularity, dimension order, etc.
    checkphysical(arrdef(:,Ne+1:end),pdims,psize);
    
%     % Get the number of B-Combiner-Boxes (potentially independent MPPT-inputs) per inverter
%     data = sortrows(data,[3 7 6]);
%     MPPTs = unique(data(:,[3,7]),'rows');
%     [~,idx,ida] = unique(MPPTs(:,1));
%     Nc = int32(diff([idx;size(MPPTs,1)+1]));
% 
%     if min(Nc)~= max(Nc)
%     warning('readpparraydef:Nccheck','readpparraydef: %d inverters have less than %d inputs', ...
%             sum(Nc~=max(Nc)),numel(Nc));
%     end
%         
%     % PROVISIONAL: ignore box1, and get MPPT-input Nº from box1 + inverter
%     %              use of the three indices requires full-implementation of NDmap
%     
%     MPPTi = 1:size(MPPTs,1);
%     % stringindex = unique(data(:,[5,1,3,2]),'rows');
%     data(:,3) = MPPTi(ida); % replace inverter column
% 
%     ArrDef = NDmap(data(:,[3,2,1]),data(:,[5,4]));

end

function s = dimsize(smth)
% Return number of dimensions of edims/pdims, whether they're a number or string
   if isnumeric(smth),s = smth; else, s = numel(smth); end
end

function [edims, edef] = checkedims(edims,edef)    
% Assume/get connection scheme, and perform checks on the connection structure

    Ne = size(edef,2);
    assert(dimsize(edims)==Ne,'checkarrdef:edimsize','Dimension mismatch');
    if isnumeric(edims), edims = ''; end
    
    if Ne == 0,return; end
    
    % Electrical dimension types:**
    etypes = {  'x: absolute index'; 
                'c: cell';
                's: string (DC series)';
                'm: module';
                'b: box (DC parallel)';
                'e: MPPT input';
                'i: inverter';
                't: Trafo';
                'a: AC busbar' };
    
    elabels = char(etypes);
    etypes = elabels(:,1)'; %'xcsmbeita'
    elabels = cellstr(elabels(:,4:end)); 

    if isempty(edims)
        assuming = true;
        switch Ne
            case 1, edims = 'x';     % absolute index
            case 2, edims = 'cs';    % cell,str
            case 3, edims = 'msi';   % mod,str,inv
            case 4, edims = 'msbi';  % mod,str,box,inv
            case 5, edims = 'msbei'; % mods,str,box,mppt,inv
            otherwise
                % Prompt the user for edims string
                if Ne < 1, error('checkarrdef:NeCase','Don''t know what to do without electrical coordinates'); end

                prompt = sprintf(['Sorry, I''m still not very good at reading headers. ' ...
                'Please compose a %d-character string out of the following code letters, ' ...
                'to let me know what the electrical array-definition means:\n\n'],Ne);
                for j = 1:numel(etypes)
                    prompt = sprintf('%s            %c: %s\n',prompt,etypes(j),elabels{j});
                end

                default = ['ms' repmat('b',1,Ne-5) 'eit'];
                edims =char(inputdlg(prompt,'Electrical Coordinates',1,{default}));

                % Check user answer: all string letters must belong to recognized types
                if isempty(edims) || numel(edims)~= Ne
                    error('checkarrdef:promptne','Wrong number of electrical dimensions');
                end
                if ~isempty(setdiff(edims,etypes))
                    error('checkarrdef:promtelts','Unrecognized electrical dimension code-letters');
                end
                
                assuming = false; % came from user
        end
    else
        assuming = false;
    end
    
    if assuming, warning('checkarrdef:assumeelec','Assuming electrical coordinates as ''%s''',edims); end
        
    % PROVISIONAL: Absolute indices should be handled as a special case, now they're just dumped
    if edims(1) == 'x'
        warning('checkarrdef:dumpx','I still don''t know what to do with indices. Dumping them...')
        [edims, edef] = checkedims(Ne-1,edef(:,2:end));  
        return; 
    else
        % Some form of pv-surface must exist at the lowest level
        if ~any(edims(1)==['m','c'])
            error('checkarrdef:PVsurf','Module or cell expected at the lowest dimension');
        end
    end
            
    % Check what can and cannot be up- and downstream of each dimension xcsmbeit
    for j = 1:Ne
        switch edims(j)
            case 'x', checkbeforeandafter(j,edims,etypes,'','~');  % must be the first
            case 'c'
                % FUTURE: cell structure allowed inside modules
                % checkbeforeandafter(j,edims,etypes,'x','~xc');  % only after x, once
                % checkbeforeandafter(j,edims,etypes,'x','sb',true);
                
                % PROVISIONAL: only meaningfull as [x]cs
                checkbeforeandafter(j,edims,etypes,'x','s',true);
                
            % FUTURE: allow series (s) and parallel (b) arrays anywhere in DC (after xc)
            % case 's', checkbeforeandafter(j,edims,etypes,'csmb','~xc'); % allow nested arrays
            % case 'b', checkbeforeandafter(j,edims,etypes,'csmb','~xc');

            % PROVISIONAL: allow series (s) only right after c/m, and boxes (b) after s/b
            case 's', checkbeforeandafter(j,edims,etypes,'cm','bei',true); 
            case 'b', checkbeforeandafter(j,edims,etypes,'sb','bei',true); % allow nested boxes

            case 'm' 
                % FUTURE: module can come after cell structure, only once
                % checkbeforeandafter(j,edims,etypes,'xcsb','~xcm'); 
                % checkbeforeandafter(j,edims,etypes,'xcsb','sbi',true);
                
                % PROVISIONAL: no support for internal cell-structure, or parallel modules
                checkbeforeandafter(j,edims,etypes,'x','s',true); 
                
            case 'e' % MPPT input must be followed by inverter
                checkbeforeandafter(j,edims,etypes,'xcsmb','ita');
                checkbeforeandafter(j,edims,etypes,'sb','i',true); 
            case 'i'
                checkbeforeandafter(j,edims,etypes,'~ita','ta'); % only once, followed by t/a
                checkbeforeandafter(j,edims,etypes,'sbme','ta',true);
            case 't', checkbeforeandafter(j,edims,etypes,'ita','ta',true); % nested trafos and arrays allowed
            case 'a', checkbeforeandafter(j,edims,etypes,'ita','ta',true);
        end
    end
end

function [pdims,pdef] = checkpdims(pdims,pdef)
% Assume/get connection scheme, and perform checks on the connection structure

    Np = size(pdef,2);
    assert(dimsize(pdims)==Np,'checkarrdef:pdimsize','Dimension mismatch');
    if isnumeric(pdims), pdims = ''; end

    if Np == 0,return; end
    
    % Physical dimension types:
    ptypes = {  'x: absolute index';
                'p: module';
                't: mount / tracker';
                'b: block' };
    
    plabels = char(ptypes);
    ptypes = plabels(:,1)'; %'xptb'
    plabels = cellstr(plabels(:,4:end));
    
    if isempty(pdims)
        assuming = true;
        switch Np
            case 0  % don't do anything (allow pure-electrical definitions)
            case 1, pdims = 'x';     % absolute index 
            case 2, pdims = 'pt';    % module,tracker
            % that's pretty much what we know...
            otherwise
                prompt = sprintf(['Yea, maybe I told you I''m not very good at reading headers. ' ...
                'Please compose a %d-character string out of the following code letters, ' ...
                'to let me know what the physical array-definition means:\n\n'],Np);
                for j = 1:numel(ptypes)
                    prompt = sprintf('%s            %c: %s\n',prompt,ptypes(j),plabels{j});
                end
                default = ['pt' repmat('b',1,Np-2)];
                pdims =char(inputdlg(prompt,'Physical Coordinates',1,{default}));

                % all string letters must belong to recognized types
                if isempty(pdims) || numel(pdims)~= Np
                    error('checkarrdef:promptnp','Wrong number of physical dimensions');
                end
                if ~isempty(setdiff(pdims,ptypes))
                    error('checkarrdef:promptplts','Unrecognized physical dimension code-letters');
                end
                assuming = false;
        end
    else
        assuming = false;
    end
    if assuming, warning('checkarrdef:assumephys','Assuming physical coordinates as ''%s''',pdims); end
    
    % PROVISIONAL: Absolute indices should be handled as a special case, now they're just dumped
    if pdims(1) == 'x'
        warning('checkarrdef:dumpx','I still don''t know what to do with indices. Dumping them...')
        [pdims, pdef] = checkpdims(Np-1,pdef(:,2:end));
        return; 
    else
        % Some form of pv-surface must exist at the lowest level
        if ~pdims(1)=='p'
            error('checkphysical:pmod','module expected at the lowest physical dimension');
        end
    end
    
    % Check what can and cannot be up- and downstream of each dimension
    for j = 1:Np
        switch pdims(j)
            case 'x', checkbeforeandafter(j,pdims,ptypes,'','~'); % must be the first
            case 'p'
                checkbeforeandafter(j,pdims,ptypes,'x','~xp'); % can come only after x, only once
                checkbeforeandafter(j,pdims,ptypes,'x','t',true); % must be followed by t
            case 't'
                checkbeforeandafter(j,pdims,ptypes,'~tb','b'); % xp before, b after
                checkbeforeandafter(j,pdims,ptypes,'p','b',true); 
            case 'b'
                checkbeforeandafter(j,pdims,ptypes,'tb','b',true); % nested blocks allowed
        end
    end
end

function flag = checkelectrical(edef,edims)
% Check uniqueness, regularity, dimension order, etc. of electrical dimensions

    [Nm,Ne] = size(edef);
    if Ne ~= numel(edims), error('checkarrdef:nedims','dimension mismatch'); end
    if Ne == 0, flag = 0; return; end

    % Check uniqueness of electrical coordinates
    if Ne > 0
        [~,idx] = unique(edef,'rows','stable');
        if numel(idx)~= Nm
            error('checkarrdef:uniqe','%d repeated electrical coordinates', Nm-numel(idx));
        end
    end
    
    % Chech that any parallel array of series of 'stuff' (cells or modules), has the same number 
    % of unit elements, i.e. same nominal voltage.
    for k = find(edims=='s')
        a = find(edims(1:k)=='m',1); % not-empty a < k represents a string of modules
        if ~isempty(a)
            % Get the dimension-index of the limiting parallel-array, with priority e > i > b, 
            % i.e. look for an independent mppt, if there isn't one get an inverter, then the 
            % highest-hierarchy connection box:
            b = bsxfun(@eq,edims(k+1:end)','bie'); % (Ne-k)·3 matrix with rows edims(j-k) == 'bie'
            [b,~] = ind2sub([Ne-k,3],find(b,1,'last'));
            b = b + k;
            
            checkstrings(a,b,k)
        else
        % if it wasn't a string of modules, look for a string of cells
            a = find(edims(1:k)=='c',1);
            assert(~isempty(a),'checkarraydefinition:strofwhaat','Expecting m/c before s in edims');
            
            % Just like for module-strings, but search above for modules, then 'boxes'
            b = bsxfun(@eq,edims(k+1:end)','bm');
            [b,~] = ind2sub([Ne-k,3],find(b,1,'last'));
            b = b + k;
            
            checkstrings(a,b,k)
        end
    end
%     % Get No. of modules-per-string: for string k, in MPPT l, Nm = mps(k,l) 
%     mps = shiftdim(sum(ArrayDef.PIDX>0,1),1); 
%     % MPPTs can have different number of strings, i.e. mps(k,l) = 0, 
%     % but not strings with different number of modules, mps(k,l) = {mps(q,l) or 0}
%     problematic = any(mps > 0 & mps ~= repmat(mode(mps,1),ArrayDef.esize(2),1),1);
%     if any(problematic)
%         switch optquestdlg(sprintf(['You really want to check the Array Definition,',...
%                 'there are strings with different numbers of modules in MPPT(s): {',...
%                 repmat('%d ',1,nnz(problematic)) '}. Do you want to exclude the MPPTs ',...
%                 'from the simulation?'],find(problematic)),'ArrayDef Check');
%             case 'Yes'
%                 excludedmppts = find(problematic);
%             case 'No'
%                 excludedmppts = [];
%                 warning('GUImain:ArrayDefCheck','Fine, but don''t come crying when you find huge missmatch losses');
%             otherwise
%                 error('Wise decision, fix it and then come back');
%         end
%     else,       
%         excludedmppts = [];
%         fprintf('\t Array Def. looks ok: %d mppts· %d strings · %d modules\n',ArrayDef.esize);
%     end

    flag = 1;
    
    % PENDING: more elaborate checks
    
    function checkstrings(a,b,k)
        if ~isempty(b)
            [~,~,ig] = unique(edef(:,b:end),'rows'); % get index of independent groups
            haveissues = zeros(size(ig));
            for j = 1:max(ig)
                [~,~,is] = unique(edef(ig == j,k:b-1),'rows'); % index of strings within group
                n = sum(bsxfun(@eq,1:max(is),is),1); % number of elements with each index
                if ~all(n == n(1))
                    haveissues(ig == j) = any(bsxfun(@eq,is,find(n ~= mode(n))),2); 
                end
            end
            if any(haveissues)
                error('checkelectrical:stuffperstring',['There are different number of %ss ',...
                        'within strings of the same %s, check arraydef %s.'],stuffname(edims(a)),...
                        stuffname(edims(b)),shortliststr(find(haveissues),'row',5));
            end
        end 
        function name = stuffname(c)
            switch c
                case 'i', name = 'inverter';
                case 'e', name = 'MPPT';
                case 'm', name = 'module';
                case 'c', name = 'cell';
                case 'b', name = 'box/parallel-array';
                otherwise, name = 'stuff';
            end
        end
    end
end

function flag = checkphysical(pdef,pdims,psize)
% flag = CHECKPHYSICAL(pdef,pdims,psize)
% flag = CHECKPHYSICAL(pdef,pdims,pidx)
% Check uniqueness, regularity, dimension order, etc. of physical dimensions
% Compare indices against maximum-index Np-vector psize, or against an 1·Np cell-array of allowed
% indices pidx.

    if iscell(psize), pidx = psize; psize = [];
    else, pidx = {}; end

    [Nm,Np] = size(pdef);
    if Np ~= numel(pdims), error('checkarrdef:nedims','dimension mismatch'); end
    if Np == 0, flag = 0; return; end
    
    % Check uniqueness of physical coordinates
    if Np > 0
        [~,idx] = unique(pdef,'rows');
        if numel(idx)~= Nm            
            error('checkarrdef:uniqp','%d repeated physical coordinates', numel(idx)-Nm);
        end
    end
    
    % Compare pdef against physical layout size.
    if ~isempty(psize)
    % NOTE the use of fliplr: array definitions are expected in order children -> parent, whereas 
    % array sizes are expected as parent -> children

        if pdims(1) == 'x' && numel(psize) == Np - 1
            psize(end+1) = max(pdef(:,1));  % add absolute-index-'size', if missing
        end  

        assert(numel(pdims) == numel(psize), 'checkarrdef:pdimsvspsize',...
            'Physical array-def. dimensions do not match provided layout-size');
        assert(all(max(pdef,[],1) <= fliplr(psize)),'checkarrdef:pdefvspsize',...
            'Physical array-def. indices exceed provided layout-size');
    end
    
    % Compare pdef against extensive list of allowed physical indices
    if ~isempty(pidx)      
        
        if pdims(1) == 'x' && numel(pidx) == Np - 1
            pidx = [{pdef(:,1)},pidx];  % add absolute-index-list if missing
        end  
        
        assert(numel(pdims) == numel(pidx), 'checkarrdef:pdimsvspidx',...
            'Physical array-def. dimensions do not match provided layout-indices size');
        
        for k = 1:Np
            assert(isempty(setdiff(pdef(:,k),pidx{k})),'checkarrdef:pdefvspidx',...
                'Physical array-def. indices do not match existing layout');
        end
    end
    
    % PENDING: more elaborate checks
    flag = 1;

end

function checkbeforeandafter(j,dims,types, beforeok,afterok,immediate)
% Take a string 'dims', in which each letter is a subset of string 'types', and check that all
% characters before dims(j) are included in the string 'beforeok', and that all characters
% after dims(j) are included in the string 'afterok'.
%
% If beforeok, or afterok start with '~...', the complement of 'types', and whatever characters 
%   follow the ~ will be used. (e.i. '~a' = all types but a, '~' = all types).
%
% If immediate flag is true (default false), only the immediate neighbors are checked.

    if nargin < 6, immediate = false; end

    if isempty(beforeok) && j ~= 1
        error('checkarrdef:firstdim','<%c> layer is expected to always be the first',dims(j));
    else
        if beforeok(1)== '~', beforeok = setdiff(types,beforeok(2:end)); end % ~ complement
        
        if immediate
            if j <= 1, idx = []; else, idx = j-1; end
        else
            idx = 1:j-1;
        end

        unwanted = setdiff(dims(idx),beforeok);
        if ~isempty(unwanted)
            if immediate, tag = 'immediately '; else, tag = ''; end
            error('checkarrdef:dimbefore','<%s> layer(s) are not expected %sbelow <%c>',unwanted,tag,dims(j));
        end
    end

    if isempty(afterok) && j ~= numel(dims)
        error('checkarrdef:lastdim','<%c> layer is expected to always be the last',dims(j));
    else
        if afterok(1)== '~', afterok = setdiff(types,afterok(2:end)); end % ~ complement

        if immediate
            if j >= numel(dims), idx = []; else, idx = j+1; end
        else
            idx = j+1:numel(dims);
        end
        
        unwanted = setdiff(dims(idx),afterok);
        if ~isempty(unwanted)
            if immediate, tag = 'immediately '; else, tag = ''; end
            error('checkarrdef:dimafter','<%s> layer(s) are not expected %sabove <%c>',unwanted,tag,dims(j));
        end
    end
end
		
