classdef NDmap
% NDMAP Class
%
% Index object class that represents a bijective relationship between two discrete vector spaces: 
% a domain E = [e,f,g,...] and a codomain P = [p,q,...]. The class was designed to represent 
% the relation between two quasi-regular trees / nested structures, i.e. discrete vector spaces in
% which each dimension is the 'parent' of the next: X(j,:) = [... grandparent, parent, child]. 
% This has important implications on the way indexing is handled (see EIDX, PIDX propperties, and
% HELP for NDmap.NDmap, NDmap.esubs, and NDmap.eidx).
%
% EXAMPLE: the connection scheme of a PV-plant can be represented by an NDMAP between a set of
%   'electrical coordinates' E = [inverter, string, module] and a set of 'physical coordinates'
%   P = [mount, position]
%
% PROPERTIES
%   esize - size vector for E space (e.g. [Ni,Ns,Nm] for Ni inverters, Ns strings, Nm modules)
%   psize - size vector for P space (e.g. [Nt,Np] for Nt mounts, Np module-positions)
%   EIDX - a fliplr(psize)-size lookup table such that EIDX(...,r,q,p) returns an index in E for   
%       the point p,q,r... in P, namely EIDX(...,r,q,p) = sub2ind(fliplr(esize),...g,f,e) so that 
%       calling ind2sub(fliplr(esize),EIDX(...,r,q,p)) returns the coordinates ...g,f,e.
%
%       NOTE: The reason to use all coordinates backwards, and fliplr(esize/psize) instead of 
%       simply esize, is so that linear indices inside EIDX are ordered [and returned] sorted as
%       parent > child. Using flipped coordinates and look-up tables ensures that a statement of 
%       the form ind2sub(fliplr(esize),EIDX(a:b,:)), intended to return the E-coordinates of all 
%       children of parents a:b, will return first all the children of a, followed by those of a+1,
%       a+2, .. up to b. For an intuitive example of the difference compare the output of:
%
%           [r{1:2}] = ind2sub([2,3],(1:6)'); [r{:}]
%           [r{[2,1]}] = ind2sub([3,2],(1:6)'); [r{:}]
%
%   PIDX - just the opposite of EIDX. 
%   Ne - number of dimensions of E-space, i.e. numel(esize)
%   Np - number of dimensions of P-space, i.e. numel(psize)
%
% METHODS
%   NDmap(elist,plist,esize,psize) - class constructor
%   esubs(obj,p,q,..) - returns the E-coordinates corresponding to the P-point p,q,r,...
%   psubs(obj,e,f,..) - opposite to psubs
%   eidx(obj,p,q,..) - returns the obj.esize-linear-index(es) for esubs(obj,p,q,..)
%   pidx(obj,e,f,..) - inverse of eidx.
%   eidxT(obj,p,q,...) - returns fliplr(obj.esize)-indices given P-coordinates.
%   pidxT(obj,p,q,...) - inverse of eidxT.
%   map2lists(obj) - Inverse of class constructor
%   mergemap(B,A) - create an NDmap such that C.psubs(args) = B.psubs(A.pidx(args)) 
%   compose(A,B,..) - create an NDmap of unique combinations of the domain/codomain of A, B,...
%   LabelMap(obj) - Type cast to subclass LabelMap
%
% METHODS (static)
%   NDmap.idx2sub(s,IND,ND) -  wrapper for ind2sub to handle vectors and collect output in array
%   NDmap.sub2idx(s,varargin) -  currently equivalent to sub2ind.
%
% See also: NDMAP.NDMAP, NDMAP.ESUBS, NDMAP.EIDX, NDMAP.EIDXT, NDMAP.MAP2LISTS, NDMAP.MERGEMAP, 
%   NDMAP.COMPOSE, NDMAP.IDX2SUB, LABELMAP

	properties (GetAccess = public, SetAccess = protected)
        esize  % size vector for E space (e.g. [Ni,Ns,Nm] for Ni inverters, Ns strings, Nm modules)
		psize  % size vector for P space (e.g. [Nt,Np] for Nt mounts, Np module-positions)
        Ne     % number of dimensions of E-space, i.e. numel(esize)
        Np     % number of dimensions of P-space, i.e. numel(psize)
		EIDX   % fliplr(psize)-size lookup table such that EIDX(...,r,q,p) = sub2ind(fliplr(esize),...g,f,e)
		PIDX   % fliplr()-size lookup table such that PIDX(...,g,f,e) = sub2ind(fliplr(psize),...r,q,p)
    end
    properties (Constant,Hidden = true)
        PRECISION = 'uint32'; % support for 4.3e9 - element arrays
    end
	methods
		function obj = NDmap(elist,plist,esize,psize)
        % obj = NDMAP(elist,plist,esize,psize) take the M·Ne and M·Np index arrays elist and plist,
        % which link M points in P to matching M points in E. i.e. 
        %       elist(j,:) = [p,q,r,...] <-> plist(j,:) = [e,f,g,...] for j = 1:M
        % If esize and psize are omitted, they will be set to max(elist,[],1) and max(plist,[],1)
        % NDmap(plist) - works as a shortcut for NDmap((1:size(plist,1))',plist)
        % NDmap() - creates an empty NDmap object (esize = psize = [0,0])
        %
        % NOTE: If the E/P spaces are used to represent trees / nested structures, the class is
        % designed considering that all indices are ordered [..., parent, child] so that function 
        % calls of the form f(parent,:) [where f stands for eidx, pidx, esubs, psubs] will return 
        % references to all the children of 'parent'.
        %
        % See also NDMAP, NDMAP.EIDX, NDMAP.ESUBS

            switch nargin
            case 0
                obj.esize = [0,0];
                obj.psize = [0,0];
                obj.Ne = 0;
                obj.Np = 0;
                obj.PIDX = zeros(0,0,NDmap.PRECISION);
                obj.EIDX = zeros(0,0,NDmap.PRECISION);
                return
            case 1
                obj = NDmap((1:size(elist,1))',elist); 
                return
            end
            
            % Allocate minimu size for E
            if nargin < 3 || isempty(esize), esize = max(elist,[],1); end 
            assert(isrow(esize),'NDmap:expecting horizontal vector esize');

            % Handle vectors specially (use two-valued sizes e.g. [n,1], for ne = 1)
            % NOTE: if matrix-style indices are used, i.e. size(elist,2) == 1, E-space will be 
            % treated as a matrix, regardles of any(esize == 1)
            if isscalar(esize) || (size(elist,2) == 1 && numel(esize) == 2 && any(esize == 1))
                ne = 1;
                if isscalar(esize), esize = [esize,1]; end  % fix for vectors
                assert(max(elist(:))<= max(esize),'NDmap:index(es) in elist out of bounds');
            else
                ne = numel(esize);
                assert(size(elist,2) == ne,'NDmap:unmatching elist and dimensions');
                assert(nnz(bsxfun(@gt,elist,esize))==0,'NDmap:index(es) in elist out of bounds');
            end

            % Do the same for P
            if nargin < 4 || isempty(psize), psize = max(plist,[],1); end
            assert(isrow(psize),'NDmap:expecting horizontal vector psize');

            if isscalar(psize) || (size(plist,2) == 1 && numel(psize) == 2 && any(psize == 1))
                np = 1;
                if isscalar(psize), psize = [psize,1]; end
                assert(max(plist(:))<= max(psize),'NDmap:index(es) in plist out of bounds');
            else
                np = numel(psize);
                assert(size(plist,2) == numel(psize),'NDmap:unmatching plist and dimensions');
                assert(nnz(bsxfun(@gt,plist,psize))==0,'NDmap:index(es) in plist out of bounds');
            end

            % Verify that there's a matching number of points in E and P
            assert(size(elist,1)==size(plist,1),'NDmap:unmatching esize and psize');
            
            % Stick to connected points (all E- and P-coordinates > 0)
            notconn = ~all([elist,plist]>0,2);
            if any(notconn)
                warning('NDmap:notconn','NDmap:%d disconnected rows ignored',nnz(notconn));
            end

            % Remove redundant rows
            idx = unique([elist,plist],'rows');
            if numel(idx) < size(elist,1)
                warning('NDmap:redundant','NDmap:%d redundant rows ignored',size(elist,1)-numel(idx));
                elist = elist(idx); 
                plist = plist(idx);
            end
            
            % Check that E- and P-coordinates are independently unique 
            idx = sort([repeatedrows(elist);repeatedrows(plist)]);
            if ~isempty(idx)
                if numel(idx)>20
                    idx = idx(1:20);
                    trimmed = true;
                else
                    trimmed = false;
                end
                msg = sprintf(['\t E:' repmat('%d ',1,ne) '<-> P:' repmat('%d ',1,np) '\n'], [elist(idx,:), plist(idx,:)]');
                msg = [sprintf(['Repeated E- and/or P-coordinates\n'...
                                'Verify the following rows in elist / plist:\n']), msg];
                if trimmed, msg = [msg, sprintf('\t...\n')]; end
                error('NDmap:graphuniqueness',msg);
            end
            
            % Create object
            obj.esize = double(esize);
            obj.psize = double(psize);
            obj.Ne = ne;
            obj.Np = np;
            
            % Get linear indices of E- and P-points
            args = fliplr(num2cell(elist,1));
			ee = NDmap.sub2idx(fliplr(obj.esize),args{:});
            args = fliplr(num2cell(plist,1));
			pp = NDmap.sub2idx(fliplr(obj.psize),args{:});
			
            % And cross-reference the two arrays EIDX / PIDX
			obj.PIDX = zeros(fliplr(obj.esize),NDmap.PRECISION);
			obj.EIDX = zeros(fliplr(obj.psize),NDmap.PRECISION);
			obj.PIDX(ee) = pp;
			obj.EIDX(pp) = ee;
            
            function idx = repeatedrows(A)
            % Returns the indices of any rows of A which are repeated

                n = size(A,1);
                [B,ndx] = sortrows(A); % sort rows of A, remeber the original order in ndx

                d = all(B(1:n-1,:)==B(2:n,:),2); % check each row with the next, see if they're equal
                d = [d;false]|[false;d];         % if d(a) = d(a+1), make both d(a:a+1) = true
                d(ndx) = d;                      % un-sort d (to the original row order)
                idx = find(d);                   % return indices
            end
        end
        
        function B = invert(A)
            B = NDmap();
            B.EIDX = A.PIDX;
            B.PIDX = A.EIDX;
            B.esize = A.psize;
            B.psize = A.esize;
            B.Ne = A.Np;
            B.Np = A.Ne;
            % B = NDmap(plist,elist,A.psize,A.esize);
        end
		
		function ee = eidx(obj,varargin)
        % ee = EIDX(obj,p,q,r..) Returns the linear index(es) of an array of size obj.esize, 
		% corresponding to the P-coordinates p,q,r,... in varargin. Subscripting with colon (:) 
        % works similar(*) to MATLAB's matrix indexing, such that, e.g. PIDX(obj,j,:) returns the  
        % indices for all points with p = j.
        % Output is an [n,1] vector of indices such that ind2sub(obj.esize,ee(j)) for j = 1:n
        % returns the E-coordinates of point j.
        %
        % CAUTION (*): internally, indexing is performed over an array of size fliplr(obj.esize), 
        % so that the resulting indices are sorted by dimensions 1:n as opposed to MATLAB's usual
        % matrix indexing (sorted n:1). 'Incomplete-indexing' (e.g. obj.eidx(1,2) when obj.Np > 2)
        % is not yet fully supported, the only exception is an ending semicolon, e.g. eidx(1,:)
        % Use a full-set of indices (e.g. obj.eidx(1,2,:,1) for Np = 4) whenever possible.
        % See help for esubs() for further detail.
        %
        % See also: NDMAP, ESUBS, PIDX
			
            % handle incomplete indexing
            if nargin - 1 > 1 && nargin - 1 < obj.Np
                if ischar(varargin{end}) && varargin{end}==':'
                    varargin(end+1:obj.Np) = {':'};
                else
                    error('NDmap:eidx','Incomplete indexing is not supported');
                end
            end
			ee = nonzeros(obj.EIDX(varargin{end:-1:1}));            
        end
			
		function pp = pidx(obj,varargin)
        % pp = PIDX(obj,e,f,g..) Analogous to eidx, except that it returns P-indices from 
        % E sub-scripts e,f,g..
        % See also: NDMAP, EIDX, PSUBS

            % handle incomplete indexing
            if nargin - 1 > 1 && nargin - 1 < obj.Ne
                if ischar(varargin{end}) && varargin{end}==':'
                    varargin(end+1:obj.Ne) = {':'};
                else
                    error('NDmap:pidx','Incomplete indexing is not supported');
                end
            end
			pp = nonzeros(obj.PIDX(varargin{end:-1:1}));
        end
        
		function elist = esubs(obj,varargin)
        % s = ESUBS(obj,p,q,r..) returns the E-coordinates corresponding to the P-point p,q,r,... 
        % in varargin. Subscripting with colon (:) works similar(*) to MATLAB's matrix indexing,  
        % such that, e.g. ESUBS(obj,j,:) returns the subscripts for all points with p = j.
        % Output is an m·Ne array representing m points in Ne-dimensional space E.
        %
        % CAUTION (*): internally indexing is performed over an array of size fliplr(obj.esize), 
        % so that the resulting indices are sorted by dimensions 1:n as opposed to MATLAB's usual
        % matrix indexing (n:1). The intention is that a call of the form obj.esubs(p,q) 
        % [where p is a vector of 'parents' and q is a vector of 'children'] will return subscripts
        % ordered by parent and then by children. Note that this is NOT how MATLAB would usually 
        % return the linear indices for a matrix (ordered by column, then by row... or in the 
        % general case by dimension end:-1:1). For an intuitive example of the difference compare 
        % the output of:
        %
        %           [r{1:2}] = ind2sub([2,3],(1:6)'); [r{:}]
        %           [r{[2,1]}] = ind2sub([3,2],(1:6)'); [r{:}]
        %
        % 'Incomplete-indexing' (e.g. obj.esubs(1,2) when obj.Np > 2) is not yet fully supported, 
        % the only exception is an ending semicolon, e.g. esubs(1,:).
        % Use a full-set of indices (e.g. obj.esubs(1,2,:,1) for Np = 4) whenever possible.
        % 
        % See also: NDMAP, EIDX, PSUBS
               
            ee = obj.eidx(varargin{:});
            elist = fliplr(NDmap.idx2sub(fliplr(obj.esize),ee,obj.Ne));
        end
        
		function plist = psubs(obj,varargin)
        % s = PSUBS(obj,e,f,..) - analogous to esubs(), except that it returns P-subscripts given
        % E-coordinates. Output is an m·Np array representing m points in Np-dimensional space P.
        % 
        % See also: NDMAP, ESUBS, PIDX
        
            pp = obj.pidx(varargin{:});
            plist = fliplr(NDmap.idx2sub(fliplr(obj.psize),pp,obj.Np));
        end
        
        function eeT = eidxT(obj,varargin)
        % eeT = EIDXT(obj,p,q..) Returns E-space indices given P-coordinates, just like eidx 
        % except that these are refered to a an esize array instead of fliplr(esize). In other 
        % words, EIDXT returns the actual linear-index corresponding to the coordinates of 
        % esubs(obj,varargin).
        %
        % See also: NDMAP, ESUBS, EIDX, PIDXT

            E = num2cell(esubs(obj,varargin{:}),1);
            eeT = NDmap.sub2idx(obj.esize,E{:});
        end
		
        function ppT = pidxT(obj,varargin)
        % ppT = PIDXT(obj,e,f..) Returns the linear indices for psubs(obj,varargin) given
        % given E-coordinates e,f,..(See HELP for EIDXT).
        %
        % See also: NDMAP, PSUBS, EIDXT, PIDX

            P = num2cell(psubs(obj,varargin{:}),1);
            ppT = NDmap.sub2idx(obj.psize,P{:});
        end
        
        function [elist,plist] = map2lists(obj)
        % [elist,plist] = MAP2LISTS(obj) 
        % Inverse of class constructor: generates the Ne- and Np-column arrays of indices that 
        % would be used to create obj - see comments on NDmap/NDmap
        %
        % See also: NDMAP, NDMAP.NDMAP
        
            conn = nonzeros(obj.EIDX(:));
            elist = fliplr(NDmap.idx2sub(fliplr(obj.esize),conn,obj.Ne));            
            plist = fliplr(NDmap.idx2sub(fliplr(obj.psize),obj.PIDX(conn),obj.Np)); 

            out = sortrows([elist,plist]);
            idx = find(all(out,2));
            elist = out(idx,1:obj.Ne);
            plist = out(idx,obj.Ne+1:obj.Ne+obj.Np);
        end
        
        function T = deftable(obj,varargin)
        % T = DEFTABLE(obj) - similar to NDMAP.MAP2LISTS (inverse of constructor), but returning 
        %   a single (flipped) table T = fliplr([psubs,esubs]).
        %   The function is provided for consistency with convention [m s i p t] to describe
        %   array-interconnection in PV-systems.
        %
        % T = DEFTABLE(obj,e,f,..) - returns a 'filtered' table, considering only the elements
        %   with electrical coordinates e,f,.... as passed to PSUBS(obj,e,f,..)
        %
        % See also: NDMAP, NDMAP.NDMAP

            if nargin == 1, varargin = {':'}; end
            pp = obj.pidx(varargin{:}); % fliplr(psize)-linear indexes
            T(:,obj.Ne+(1:obj.Np)) = NDmap.idx2sub(fliplr(obj.psize),pp,obj.Np);            
            T(:,1:obj.Ne) = NDmap.idx2sub(fliplr(obj.esize),obj.EIDX(pp),obj.Np); 
        end
        
        function C = mergemap(B,A)
        % C = MERGEMAP(B,A)
        % Take NDmap objects A and B such that all(A.psize <= B.esize) and compose a third NDmap C
        % such that C.psubs(args) = B.psubs(A.pidx(args)) and C.esubs(args) = A.esubs(B.eidx(args))
        %
        % See also: NDMAP.COMPOSE, NDMAP
            
            assert(isa(A,'NDmap') && isa(B,'NDmap'),'NDmap_merge:args','Expecting two NDmap objects');
            assert(all(A.psize <= B.esize),'mergemap:size','Incompatible NDmaps');
            
            [Ae,Ap] = map2lists(A);
            [Be,Bp] = map2lists(B);
            
            % Get elements that are connected in both A.p and B.e:
            [connA,locB] = ismember(Ap,Be,'rows');
            
            % Create new map merging connected E-points from A to P-points in B
            C = NDmap(Ae(connA,:),Bp(nonzeros(locB),:),A.esize,B.psize);
        end
        
        function M = compose(varargin)
        % M = COMPOSE(A,B,C,...) Take NDmap objects A, B, C... and compose an NDmap M whose 
        % domain and codomain are all unique combinations of the domain and codomain of A, B, C...
        % M.esize will be a concatenation of non-singleton [A.esize, B.esize, C.esize...] and 
        % similarly for M.psize.
        %
        % See also: NDMAP.MERGEMAP, NDMAP
        
            if nargin < 1, M = NDmap(); return; end
            
            for j = nargin:-1:1
                X = varargin{j};
                [E{j},P{j}] = map2lists(X);     % get point lists for each NDmap
                n(j) = size(E{j},1);            % save point number
                s{1,j} = X.esize(1:X.Ne);       % save esize and psize,
                s{2,j} = X.psize(1:X.Np);       % ... use of 1:Ne removes trailing ones
            end
            
            % Generate a set of all possible combinations of points
            ii = cell(1,numel(n));
            [ii{:}] = ind2sub(n(:)',1:prod(n));
            ii = reshape([ii{:}],[],numel(n));
                    
            E = arrayfun(@(j) E{j}(ii(:,j)),1:nargin,'unif',0); E = [E{:}];
            P = arrayfun(@(j) P{j}(ii(:,j)),1:nargin,'unif',0); P = [P{:}];
            
            es = [s{1,:}]; ps = [s{2,:}];   % Get combined dimensions
            M = NDmap(E,P,es,ps);           % Create composed NDmap
        end
        
        function B = LabelMap(A)
        % B = LABELMAP(A) - Type-cast to subclass LABELMAP
        % See also: LABELMAP
        
            B = LabelMap();
            for f = {'EIDX','PIDX','esize','psize','Ne','Np'}
                B.(f{:}) = A.(f{:});
            end
        end
        
        function obj = esubmap(obj,varargin)
        % SUB = ESUBMAP(OBJ,e,f,..) - Reduce NDMAP OBJ's domain E to E(e,f,..), keeping the number
        %   of dimensions (OBJ.Ne and OBJ.Np) and the size of the codomain P (OBJ.psize) intact.
        %   Points of P no longer referenced by the new filtered E will be 'disconnected' (their
        %   indices set to zero).
        
            varargin(end+1:obj.Ne) = {':'};
            varargin = fliplr(varargin);

            % subset PIDX
            obj.PIDX = obj.PIDX(varargin{:});
            obj.esize = fliplr(arrayfun(@(d) size(obj.PIDX,d),1:obj.Ne));

            % ... and rebuild EIDX
            obj.EIDX(:) = 0;
            obj.EIDX(nonzeros(obj.PIDX)) = find(obj.PIDX(:));
        end
        
        function obj = psubmap(obj,varargin)
            varargin(end+1:obj.Np) = {':'};
            varargin = fliplr(varargin);
            obj.EIDX = obj.EIDX(varargin{:});
            obj.psize = fliplr(arrayfun(@(d) size(obj.EIDX,d),1:obj.Np));
            obj.PIDX(:) = 0;
            obj.PIDX(nonzeros(obj.EIDX)) = find(obj.EIDX(:));
        end
        
        function varargout = subsref(obj,s)
           varargout = {builtin('subsref',obj,s)}; 
        end
    end
    
    methods (Static)
        function idx = sub2idx(siz,varargin)
        % idx = SUB2IDX(siz,varargin) - currently equivalent to sub2ind
        % See also: SUB2IND
        
            %if numel(siz) == 1, siz = [siz,1]; end
            idx = sub2ind(siz,varargin{:});
        end
        function subarr = idx2sub(siz,IND,ND)
        % subarr = IDX2SUB(siz,IND,nd)
        % Wrapper for ind2sub to handle vectors and collect output in a single array.
        % Unlike ind2sub (that returns a row vector when varargin is a row vector) IDX2SUB always
        % returns a vertical vectors when siz = [1,n] or [n,1], EXCEPT if explicitly ND == 2
        % subarr: is a single array, equivalent to [varargout{:}] of ind2sub.
        %
        % See also: IND2SUB
            if nargin < 3, ND = nnz(siz > 1); end
            
            if ND <= 1 && numel(siz) == 2
                subarr = ind2sub(siz,IND);
                subarr = subarr(:);
            else
                out = cell(1,numel(siz));
                [out{:}] = ind2sub(siz,IND);
                subarr = reshape([out{:}],[],numel(siz));
            end
        end
    end
end
