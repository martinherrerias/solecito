classdef LabelMap < NDmap
% LABELMAP < NDmap
%
% Index object between an esize array of integers E and a psize array of string labels P, that is
% E(j) = [e,f,g,...] <-> P(j) = {'p','q',...} for a set of points j (domain/codomain) in each 
% space. e.g. P(j) = E(j) = [1,20] <-> ['inverter_A', 'box_1', 'string_x']
% Just like NDmap, The class is designed to represent the relation between two trees / nested 
% structures, in wich each dimension is the 'parent' of the next (see NDmap for details).
% LABELMAP contains an NDmap object, plus a cell-array of string labels for each index in each
% dimension of the P-space of the NDmap.
%
% PROPERTIES
%     Lbl - Np cell-array of P dimension labels. For each dimension j, Lbl{j} is a psize(1:j)  
%       cell-array of labels. This allows for children of the same parent to have different 
%       labels for the same index. e.g. E{[1,1;1,2]} <-> P{['p','qa';'p','qb']}
%
% (Inherited from NDmap)
%     esize - size vector for E space
%     psize - size vector for P space
%     EIDX - fliplr(psize)-size lookup table of linear indices of E-space
%     PIDX - fliplr(esize)-size lookup table of linear indices of P-space
%     Ne - number of dimensions of E-space, i.e. numel(esize)
%     Np - number of dimensions of P-space, i.e. numel(psize)
%
% METHODS
%     LabelMap(elist,plist,esize,psize) - class constructor
%     plabels(obj,e,f..) - returns the label-array corresponding to the P-point(s) p,q,r,...
%
% (Overloaded)
%     eidx(obj,p,q,..) - like NDmap.eidx, but can take p-labels as input coordinates
%     mergemap(A,B) - analogous to NDmap.mergemap, where B is a LabelMap object 
%
% (Inherited)
%     esubs(obj,p,q,..) - returns the E-coordinates corresponding to the P-point p,q,r,...
%     psubs(obj,e,f,..) - opposite to psubs
%     pidx(obj,e,f,..) - inverse of eidx.
%     eidxT(obj,p,q,...) - returns fliplr(obj.esize)-indices given P-coordinates.
%     pidxT(obj,p,q,...) - inverse of eidxT.
%     map2lists(obj) - Inverse of class constructor - returns indices, not labels
%     compose(A,B,..) - create an NDmap of unique combinations of the domain/codomain of A, B,...
% 
% METHODS (static)
%   [B,L] = uniqueidx(A,s) - From cell-array of labels A (representing a tree), generates a matrix  
%        of minimal indices B, and a [1,numel(s)] cell-array of dimension-labels L.
%
% See also: NDMAP, LABELMAP.LABELMAP, LABELMAP.PLABELS, LABELMAP.EIDX, LABELMAP.MERGEMAP, 
%   LABELMAP.UNIQUEIDX


	properties
        % See NDmap for inherited propperties
        Lbl     % Np cell-array of P dimension labels. For each dimension j, Lbl{j} is a psize(1:j)  
                % cell-array of labels.
	end
	methods
		function obj = LabelMap(varargin)
        % obj = LABELMAP(elist,plist,esize,psize)
        % take the M·Ne index array elist and the M·Np cell-array of strings plist , which link 
        % M points in E to matching M points in P. 
        % i.e. elist(j,:) = [e,f,g,...] <-> plist(j,:) = {'str_p','str_q','str_r',...}
        % If esize is omitted, it will be set to max(plist,[],1)
        % If psize is omitted, it will be set to fit the maximum number of unique element-labels 
        %   at every given dimension of P. 
        % LabelMap(plist) - works as a shortcut for LabelMap((1:size(plist,1))',plist)
        % LabelMap() - creates an empty LabelMap object (esize = psize = [0,0], Lbl = {})
        %
        % NOTE: Just like NDmap, the class is designed to work for tree-like spaces, considering t
        % that all indices are ordered [..., parent, child].
        %
        % See also NDMAP.NDMAP, LABELMAP.EIDX, LABELMAP.PLABELS

            if nargin == 0
                lbls = {};
                args = {};
            else
                if nargin == 1
                    plist = varargin{1};
                    elist = (1:size(plist,1))'; 
                else
                   elist = varargin{1};
                   plist = varargin{2};
                end
                
                if nargin < 4, psize = []; else, psize = varargin{4}; end
                
                % Make sure numbers are stored as strings - unique() won't work otherwise
                plist(:) = arrayfun(@(j) num2str(plist{j}),1:numel(plist),'unif',0);
                [pidxlist,lbls] = LabelMap.uniqueidx(plist,psize);
                args = {elist,pidxlist};
            end

            obj = obj@NDmap(args{:},varargin{3:end}); % Call to superclass constructor
            obj.Lbl = lbls;
        end
		
		function ee = eidx(obj,varargin)
        % ee = EIDX(obj,'p',q,{r}..)
        % When arguments in varargin are integers or ':', EIDX works exactly like NDmap.EIDX (*)
        % However, dimension labels can also be passed as arguments (as strings or cell-arrays),
        % in which case the labels will be searched-for in the corresponding element of obj.Lbl, 
        % and only matching indices will be returned. Examples:
        %
        %   EIDX(obj,1,:) - returns all elements in E such that P = [1,q,r...], or more strictly
        %       P = {obj.Lbl{1}{1},q,r...}.
        %       (*) NOTE that numeric indices are generated by order of appearance of their label
        %       in LabelMap(). When a specific element is required, it's a better idea to use
        %       exact labels. indices might be more efficient inside loops, where the order of
        %       evaluation doesn't really matter.
        %   EIDX(obj,'A',:,'B') - return all elements in E such that P = {'A',q,'B',...} for any q
        %   EIDX(obj,{'A','B'},:) - return all elements in E such that P = {'A',...} OR {'B',...}
        %   EIDX(obj,1,:,{'A,'B'}) - combination of integer indices and labels is possible
        %
        % NOTE: since labels can be different inside different elements of the same 'dimension',
        %   label-indexing might not produce the results you would expect from integer-indexing. 
        %   E.g. for M = LabelMap({'A','B';'C','D'}), M.EIDX('A','D') returns an empty array,
        %   because there is no label 'D' for elements {'A',..} even though 'D' is a valid label
        %   for 'C' - M.EIDX('C','D') = 2.
        %   Also, as with NDmap.eidx, 'Incomplete-indexing' is still not fully supported.
        %
        % See also: NDMAP.EIDX, LABELMAP, LABELMAP.PLABELS, NDMAP.ESUBS, NDMAP.PIDX

            % handle incomplete indexing
            if nargin - 1 > 1 && nargin - 1 < obj.Np
                if ischar(varargin{end}) && varargin{end}==':'
                    varargin(end+1:obj.Np) = {':'};
                else
                    error('LabelMap:eidx','Incomplete indexing is not supported');
                end
            end

            cellnum2str = @(x) arrayfun(@(j) num2str(x{j}),1:numel(x),'unif',0);
            
            % Combined label and conventional indexing works by checking which dimensions are
            % constrained by labels, and which by numbers/semicolon...
            
            idxargs = varargin;
            islabel = false(1,numel(varargin));
            for k = 1:numel(varargin)
                % Assume strings other than ':' are labels, and turn them into cell-arrays
                if ischar(varargin{k})
                    if ~strcmp(varargin{k},':'), varargin{k} = varargin(k); end
                end
                
                if iscell(varargin{k})
                    islabel(k) = true;
                    varargin{k}(:) = cellnum2str(varargin{k}); % the same was done in constructor
                    
                    % For a first filter, relax all labeled-dimension-constraints
                    idxargs{k} = ':';
                else
                    islabel(k) = false;
                end
            end
      
            % Consider all integer/semicolon constraints
            ee = nonzeros(obj.EIDX(idxargs{end:-1:1})); 
            
            if any(islabel)
                % ... get the array of labels for the filtered points
                pcell = obj.plabels(ee);
                filter = true(numel(ee),1);
                for k = find(islabel)
                % ... and see which comply with all labeled-constraints
                    if ~any(filter), break; end
                    dimfilter = false(nnz(filter),1);
                    for j = 1:numel(varargin{k})
                        dimfilter = dimfilter | strcmp(varargin{k}(j),pcell(filter,k));
                    end
                    filter(filter) = dimfilter;
                end
                ee = ee(filter);
            end
        end
        
		function pcell = plabels(obj,varargin)
        % pcell = PLABELS(obj,e,f,...) - Returns the P-coordinate-labels corresponding to the 
        % E-point(s) e,f,g,... Subscripting works as with NDmap.psubs, such that, 
        % e.g. PLABELS(obj,j,:) returns the subscript labels for all points with p = j.
        % as with NDmap.psubs, 'Incomplete-subscripting' is not fully supported.
        % Output is an m·Np cell-array of strings, representing m points in Np-dimensional space P.
        %
        % See also: NDMAP.PSUBS, NDMAP.EIDX, LABELMAP, LABELMAP.LABELMAP

            plist = psubs(obj,varargin{:});
            if isempty(obj.Lbl)
                pcell = reshape(cellstr(num2str(plist(:))),size(plist));
                return;
            end
            
            pcell = cell(size(plist));
            plist = num2cell(plist,1); % split in cell-array of column arguments for sub2ind
            for k = 1:size(plist,2)
                idx = sub2ind([obj.psize(1:k),1],plist{1:k});
                pcell(:,k) = obj.Lbl{k}(idx);
            end
        end
        
        function C = mergemap(B,A)
        % C = MERGEMAP(B,A) - Take NDmap/LabelMap object A and LabelMap object B such that 
        % all(A.psize <= B.esize) and compose a third NDmap C such that 
        % C.plabels(args) = B.plabels(A.pidx(args)) and C.esubs(args) = A.esubs(B.eidx(args))
        %
        % See also: NDMAP.MERGEMAP, NDMAP.COMPOSE, LABELMAP
            
            C = LabelMap(mergemap@NDmap(B,A));  % merge internal NDmaps
            C.Lbl = B.Lbl;                      % ... and add labels
        end
        
%         function varargout = subsref(obj,s)
%            varargout = {builtin('subsref',obj,s)}; 
%         end
    end
    
    methods (Static)
        function [B,L] = uniqueidx(A,s)
        % [B,L] = UNIQUEIDX(A,s)
        % From cell-array of labels A, generate a matrix of minimal indices B, and a [1,n]-
        % cell-array of dimension-labels, assuming that A represents a nested tree (in which
        % elements of the first dimensions contain sub-elements of the following).
        %
        % n: number of dimensions, or 'depth' of A, n == size(A,2) == numel(s).
        % s: maximum dimension size, e.g. s = [2,3] means 2 'parents' (different labels in 
        %    A(:,1) with max. 3 'children' (different labels in A(:,2) with the same 'parent').
        % B: matrix of minimal indices. B(:,1) will range from 1 to the number of parents,
        %    B(:,2) will range from 1 to the maximum number of children per parent, B(:,3) will
        %    go from 1 to the max. number of grand-children per children, and so on.
        % L: n cell-array of dimension labels. For each dimension j, L{j} is a s(1:j) cell-array 
        %    of the labels in A, such that they match with non-empty elements in B.
        %
        % See also: LABELMAP, LABELMAP.LABELMAP

            B = ones(size(A));
            n = size(A,2);
            L = cell(1,n);

            % List of 'parents' is straightforward
            [L{1},~,B(:,1)] = unique(A(:,1),'stable');
            
            % For all other dimensions
            for k = 2:n
                % Get group-indices for all existing combinations of higher dimensions
                % NOTE use of uniquecell, as MATLAB's unique() doesn't work with cells and 'rows'
                [~,~,ci] = uniquecell(A(:,1:k-1),'rows','stable');
                for j = 1:max(ci)
                    % Get minimal indices within each group
                    [~,~,B(ci == j,k)] = unique(A(ci == j,k),'stable');
                end
            end
            if isempty(s)
                s = max(B,[],1); % Get minimal size when missing
            else
                % Check size if provided externally
                assert(all(s >= max(B,[],1)),'LabelMap:psize','Labels cannot fit in psize array');
            end

            % Fill-in lower dimension labels recursively...
            for k = 2:n
                L{k} = assignidx(A(:,1:k),s(1:k));
            end
            
            function L = assignidx(A,s)
            % Recursive assignation of unique elements of cell-array A in a cell-array of size s
                L = cell([s,1]);
                [u,~,c] = unique(A(:,1),'stable');
                if size(A,2) == 1
                    L(1:numel(u)) = u;
                else
                    for i = 1:s(1)
                        l = assignidx(A(c == i,2:end),s(2:end));
                        L(i,:) = l(:);
                    end
                end
            end
        end     
    end
end
