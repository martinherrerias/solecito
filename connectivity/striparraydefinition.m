function [arrdf,arridx] = striparraydefinition(arraydef,edims,pdims,tidx,pidx)
% Ensure that all indices in arraydef are minimal, i.e. that they range from 1 to N, where N is the 
% number of elements actually connected.
% Store the original indices as a structure of look-up tables (NDmap-objects), such that, in the
% simplest case (exceptions below), for each dimension j, arraydef(:,j) = arridx.j.pidx(arrdf(:,j))
%
% PROVISIONAL: The function also removes any 'unhandled' dimensions (anything other than msi@pt),
% in these cases, removed dimensions will be collapsed into existing ones:
%
%          arraydef(:,m) = arridx.m.psubs(arrdf(:,m)) - modules
%   arraydef(:,[e..b,s]) = arridx.s.psubs(arrdf(:,s)) - boxes and MPPT's into strings
%   arraydef(:,[a..t,i]) = arridx.i.psubs(arrdf(:,i)) - trafos and busbars into inverters
%
% Physical dimensions can be linked to non-minimal mount and position indices (tidx, pidx), if 
% these are provided, while creating the look-up tables, it will be verified that every (p,t) is 
% included in tidx,pidx
%          arraydef(:,p) = arridx.p.psubs(arrdf(:,p)) - module positions
%    arraydef(:,[..b,t]) = arridx.t.psubs(arrdf(:,t)) - trafos and busbars into inverters
%
% NOTE: that dimensions are reversed (parent->child) in NDmaps, so that they better represent a 
%   tree structure, i.e. arridx.s.eidx(j,:) returns all strings in MPPT-input j
%
% arraydef,edims, and pdims are expected to be the output of checkarraydefinition()
% tidx,pidx should be generated during layout and mount-geometry importing, and included in Trackers

    if nargin < 2 || (size(arraydef,2)==5 && isempty(edims)), edims = 'msi'; end
    if nargin < 3 || (size(arraydef,2)==5 && isempty(pdims)), pdims = 'pt'; end
    if nargin < 4, tidx = []; elseif isrow(tidx), tidx = tidx'; end
    if nargin < 5, pidx = []; elseif isrow(pidx), pidx = pidx'; end

    % At least three electrical coordinates (module, string, inverter) plus two physical 
    % coordinates (mount, position) are required:
    if ~all(any(bsxfun(@eq,'msi',edims(:))))
       error('striparrdef:msi','Array definition requires at least [m,s,i] electrical coordinates');
    end
    if ~all(any(bsxfun(@eq,'pt',pdims(:))))
       error('striparrdef:pt','Array definition currently requires at least [p,t] physical coordinates');
    end
    
    % find key array dimensions
    dimidx.m = find(edims == 'm');
    dimidx.s = find(edims == 's');
    dimidx.i = find(edims == 'i');
    dimidx.p = find(pdims == 'p') + numel(edims);
    dimidx.t = find(pdims == 't')+ numel(edims);
    
    % Make index for modules
    arridx.m = makeidx(arraydef(:,dimidx.m));

    % Make index for strings, including connection boxes ('b') and MPPT inputs ('e')
    dimidx.be = find(any(bsxfun(@eq,'eb',edims(:)),2))';
    if ~isempty(dimidx.be)
        assert(all((dimidx.be > dimidx.s) & (dimidx.be < dimidx.i)),'striparrdef:be','b/e dimensions out of place');
        dimidx.s = [dimidx.s,dimidx.be];
    end
    arridx.s = makeidx(arraydef(:,dimidx.s));
    
    % Make index for inverters, including trafos and AC busbars
    dimidx.ta = find(any(bsxfun(@eq,'ta',edims(:)),2))';
    if ~isempty(dimidx.ta)
        assert(all(dimidx.ta > dimidx.i),'striparrdef:ta','t/a dimensions out of place');
        dimidx.i = [dimidx.i,dimidx.ta];
    end
    arridx.i = makeidx(arraydef(:,dimidx.i));
    
    % Handle physical coordinates more carefully, as they're linked to the Trackers structure
    if isempty(pidx)
        pidx = arraydef(:,dimidx.p);
    else
        assert(isempty(setdiff(arraydef(:,dimidx.p),pidx)),'striparrdef:pidx',...
                                                       'Module-position index(es) out of bounds');
    end
    arridx.p = makeidx(pidx(:));
    
    % Mount indices might (in the FUTURE) include nested-block dimensions (ptbbb...), so check the
    % set of index rows, allowing the possibility for an Ntr·m - sized tidx.
    dimidx.bb = find(pdims=='b') + numel(edims);
    if ~isempty(dimidx.bb)
        dimidx.t = [dimidx.t,dimidx.bb];
    end
    
    if isempty(tidx)
        tidx = arraydef(:,dimidx.t);
    else
        assert(isempty(setdiff(arraydef(:,dimidx.t),tidx,'rows')),'striparrdef:tidx',...
                                                              'Mount index(es) out of bounds');
    end
    arridx.t = makeidx(tidx);

    % Generate reduced array definition, replacing all indices
    arrdf(:,1) = idxlist(arridx.m,dimidx.m);
    arrdf(:,2) = idxlist(arridx.s,dimidx.s);
    arrdf(:,3) = idxlist(arridx.i,dimidx.i);
    arrdf(:,4) = idxlist(arridx.p,dimidx.p);
    arrdf(:,5) = idxlist(arridx.t,dimidx.t);
    
    function l = idxlist(obj,dd)
    % Generate vector of linear E-indices from (inverted) array of P-subs
    % FUTURE: evaluate the necessity to add this functionality to NDmap
        args = num2cell(arraydef(:,fliplr(dd)),1);
        ee = sub2ind(obj.psize,args{:});
        l = obj.EIDX(ee);
    end

    function idx = makeidx(arr)
    % Create a minimal index for unique rows of arraydef(:,edims) 
        C = unique(fliplr(arr),'rows'); 
        idx = NDmap((1:size(C,1))',C);
    end
end
