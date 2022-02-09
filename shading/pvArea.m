classdef pvArea < matlab.mixin.Copyable
% PVAREA - Defines a tree structure of polygonal areas, meant to represent photovoltaic elements  
% in a 'uniform' 2D-array, structured according to a hierarchical connection scheme.
% 'Uniform' means that all elements (children) of a PVAREA are expected to be (translated) copies 
% of a single underlying PVAREA. All elements must have the same number of sub-elements; all 
% branches must thus have the same depth; and it should thus be possible to arrange the deepest-
% elements (leaves) of a single object A over a regular (A.depth-1)-dimensional, (A.dims)-size array.

    properties (GetAccess = public, SetAccess = public, AbortSet, SetObservable)
        border      % polygon object (structure with vectors .x .y)
		elements    % array of pointers to smaller pvArea objects (children)
    end
    properties (Transient = true, GetAccess = public, SetAccess = private) %(*)
		area		% border polygon area
		dims		% obj.depth-1 vector, with maximum number of elements at each depth
    end
    properties (Dependent = true)
        depth       % longest path from root to smallest child, 1 = PVAREA with no elements
    end
    methods(Access = protected)
       function cp = copyElement(obj)
       % Called by copy(obj), recursive copy of obj, including all children elements.
		
            if ~isempty(obj.elements)
            % make a recursive copy of all elements
                ee = arrayfun(@copy,obj.elements);
            else
                ee =  pvArea.empty;
            end
            cp = pvArea(obj.border,ee);
        end 
    end
    methods
        function obj = pvArea(varargin)
        % OBJ = PVAREA()returns a dot at the origin (empty area)
        %
        % OBJ = PVAREA(S) where S is a (nested) structure with fields {border,elements} (e.g. the 
        %   result of S = STRUCT(A,..) from a PVAREA object A), attempts to convert the structure 
        %   (back) into a PVAREA object. Any packed polygons will be unpacked.
        %
        % OBJ = PVAREA(X,..) creates a PVAREA object with OBJ.border = POLYGON(X,..), and no child
        %   elements. See POLYGON for details on the interpretation of the arguments.
        %
        % OBJ = PVAREA(X,..,E1,..EN) where E1,..,EN are PVAREA elements creates a PVAREA object 
        %   with OBJ.border = POLYGON(X,..), and elements E1,..EN.
        
            if nargin == 1 && isstruct(varargin{1})
               obj = pvArea.loadobj(varargin{1});
               return;
            end
        
            ispvarea = cellfun(@(x) isa(x,'pvArea'),varargin);
            if ~any(ispvarea)
                e =  pvArea.empty; % empty array of pvArea objects
                if ~isempty(varargin)
                    varargin = varargin(~cellfun(@isempty,varargin));
                end
            else
                e = varargin{ispvarea};
                varargin = varargin(~ispvarea);
            end
        
            obj.elements = e;
            
            % see documentation of polygon() for interpretations of x,y
            if ~isempty(varargin) || isempty(e)
                if isscalar(varargin) && isa(varargin{1},'polygon')
                    obj.border = varargin{1};
                else
                    obj.border = polygon(varargin{:});
                end
            else
                obj.border = areaenvelope(obj);
            end

            obj = update(obj); % update transient properties
            
            addlistener(obj,'border','PostSet',@obj.update);
            addlistener(obj,'elements','PostSet',@obj.update);
        end

        function obj = update(obj,prop,~)
        % OBJ = LOADOBJ(OBJ) -  Parse user-accessible properties to create/update transient properties:
        %   r0, cX, fun, c0, and invf.
        %
        % OBJ = LOADOBJ(OBJ,PROP) - (where PROP is a META.PROPERTY object) is meant to be called by the 
        %   event-listeners set up by the constructor upon user-accessible variables: normalize, clip,
        %   prjfun, tol, and fiam. PROP.Name is used to run only the sections of LOADOBJ that are
        %   affected by a change of PROP.
        %
        % See documentation for POLYPROJECTOR.POLYPROJECTOR.

            if nargin < 2, prop = []; end
            changed = @(s) isempty(prop) || any(strcmp(prop.Name,s));
            
            if changed('border'), obj.area = obj.border.area(); end
            if changed('elements'), obj.dims = getdims(obj); end
        end
        
		function [s,varargout] = getdims(obj,recursive)
        % S = GETDIMS(OBJ), (re)calculate OBJ.dims, i.e. an (OBJ.depth - 1) vector with the
		% 	maximum number of elements-per-parent at each hierarchical level.
		%   As an example, S = [3,2,4] means that OBJ has 3 elements, each with 2 sub-elements 
        %   (or less), each of those with 4 sub-sub-elements (or less).
        %
		% [n1,n2,n3...] = GETDIMS(OBJ) returns scalar element numbers at each level, filling in
        %   with nulls if nargout > OBJ.depth - 1.
        %
		%   NOTE: unlike MATLAB's size(), if nargout < OBJ.depth - 1, remaining dimensions are 
        %       not returned on the last element.
		%
        % .. = GETDIMS(OBJ,RECURSIVE) - By default, GETDIMS assumes that the transient propperty
        %   OBJ.elements(j).dims is up-to-date for all elements j of OBJ, and uses these values
        %   to return the parent's dimensions. The object listeners set-up by PVAREA upon object
        %   construction should ensure that this is always the case. For debugging / OCD users
        %   setting RECURSIVE to true makes sure to traverse the complete hierarchy.
        
            if nargin < 2, recursive = false; end
            n = numel(obj.elements);
        
			% return 0 for empty objects, and 1 for a single, indivisible pvArea
			if n == 0, s = []; return; end

            if recursive
            % element sizes, recalculating at each step (slow!)
                e = arrayfun(@(x) getdims(x,true),obj.elements,'unif',0);
            else
            % element sizes, trusting that they hold the correct obj.dims tag
                e = {obj.elements.dims};
            end
            d = cellfun(@numel,e);                          % element depth 
            D = max(d);
            if D == 0, s = n; return; end   
            
            short = d < D;
            if any(short)  
                warning('pvArea:nedims','Non-uniform array, returning max. size');
                e(short) = cellfun(@(x) [x,zeros(1,D-numel(d))],e(short),'unif',0);
            end
            
            % and return n followed by the maximum of each column (i.e. each depth level)
            s = [n,max(cat(1,e{:}),[],1)];
            
            if nargout > 1
                varargout = num2cell(s(2:end));
                varargout{end+1:nargout} = 0;
                s = s(1);
            end
		end
		        
		function d = get.depth(obj), d = numel(obj.dims) + 1; end
        
        function B = splitArea(A,m,n)
        % Split a Rectangular(!) pvArea 'A' into m·n smaller elements.
        % if A.depth > 1 (i.e. A is made up of smaller elements) then these are split into elements
        % one depth level below.
        
            B = A; % shallow copy of A
			
			% if A is made of smaller elements, split recursively
			if A.depth > 1
                B.elements = arrayfun(@(x) splitArea(x,m,n),B.elements);
			else
				% split this (single) pvArea
				xlim = [min(B.border.x),max(B.border.x)];
				ylim = [min(B.border.y),max(B.border.y)];
                if ~all(B.border.x==xlim(1)|B.border.x==xlim(2))||...
				   ~all(B.border.y==ylim(1)|B.border.y==ylim(2))
					error('pvarea:split','Cannot split non-rectangular areas');
                end
                
                % create a small rectangular element (hard copy of new element)
				w = diff(xlim)/n;
				h = diff(ylim)/m;
                b = pvArea(w,h);

                % make each element of b a new pvArea, translated by x0, y0
                [x0,y0] = meshgrid(xlim(1)-w/2+(1:n)*w,ylim(1)-h/2+(1:m)*h);
                B.elements = arrayfun(@(x,y) translate(copy(b),x,y),x0,y0);
			end
		end

        function A = areaArray(a,xc,yc) 
		% areaArray(a,xc,yc) replicates pvArea 'a' along positions xc,yc and saves the result
		%		as the elements of pvArea A. 
        % areaarray(a,[m n]) creates a grid of m rows, n columns; with minium spacing between
		%		rows and columns as determined by limits of a.x and a.y.
		%		The resulting array is centered around 0,0.
            
            narginchk(2,3);
            if nargin == 2
                if numel(xc)~=2
                    error('areaarray:input','use either a,xx,yy or a,[m,n] as input parameters!');
                else
                    m = xc(1); 
                    n = xc(2);
                    w = max(a.border.x)-min(a.border.x);
                    h = max(a.border.y)-min(a.border.y);
                    [xc,yc] = meshgrid((-0.5*(n-1):0.5*(n-1))*w,(-0.5*(m-1):0.5*(m-1))*h);
                end
            else
               assert(numel(xc)==numel(yc),'expecting equal-sized xc, yc'); 
            end
            
            A = pvArea();			   % initialize empty A

            % replicate element 'a' with all offsets xc, yc
            A.elements = arrayfun(@(x,y) translate(copy(a),x,y),xc,yc);
            A.border = areaenvelope(A);
        end

        function obj = translate(obj,x,y)
		% OBJ = TRANSLATE(OBJ,X,Y) - Translate obj and all its elements by X,Y
        
            narginchk(2,3);
            if nargin < 3
                assert(numel(x) == 2,'Expecing 2-vector or X, Y argument pair');
                y = x(2); x = x(1);
            else
                assert(isscalar(x) && isscalar(y),'Expecing 2-vector or X, Y argument pair');
            end
		
			% apply to root object
			obj.border = polytranslate(obj.border,[x,y]);
			if obj.depth > 1
				% apply recursively to each element
                obj.elements = arrayfun(@(e) translate(e,x,y),obj.elements);
			end
        end
        
        function obj = rotate(obj,th)
		% Rotate obj.border and all its elements by an angle th (degrees)
		
			% apply to root object
			obj.border = polyrotate(obj.border,th);
			if obj.depth > 1
				% apply recursively to each element
                obj.elements = arrayfun(@(e) rotate(e,th),obj.elements);
			end
        end
		
        function B = flattenpvarea(A)
        % B = FLATTENPVAREA(A) - Returns a pvArea of depth 2 where B.elements is a list of the 
        % deepest sub-elements (leaves) of the original pvArea object A.
			
            switch numel(A)
            case 0, B = pvArea.empty();
            case 1
                if isempty(A.elements), B = copy(A); return; end
                
                E = A.elements;
                d = E(1).depth;
                assert(all([E.depth] == d),'FLATTENPVAREA requires uniform depth');
                for j = 1:d-1
                   E = [E.elements];
                end
                B = pvArea(A.border,E);
            otherwise
                error('FLATTENPVAREA requires scalar object');
            end
		end
		
		function newobj = reshapearea(obj,arraydef)
        % NEWOBJ = RESHAPEAREA(OBJ, ARRAYDEF)
		% Take the N elements of pvArea object obj and re-arrange them according to N·k index array 
		% arraydef, such that:
		%
		%	newobj.elements(a(j,k)).elements(a(j,k-1))...elements(a(j,1)) = obj.elements(j) [or obj(j)]
		%	where a(j,l) is short for arraydef(j,l)
		%
		% The result is a pvArea of depth k-1+obj.depth, inside which the original N elements of obj are
		% arranged as sub-elements
		%
		% Note that if numel(arraydef) == N (i.e. single column), it will just be used as a new sorting
		% order for the elements of obj.
		%
		% If an array of N pvArea objects is provided as obj, they will first be grouped in as elements of
		% a parent object.
		%
		% arraydef can also be a 3-dimensional array, of size nr·nc·k such that nr·nc == N
		%			the array will just be reshaped to work as usual.

            if numel(obj) > 1
				newobj = reshapearea(pvArea(obj(:)),arraydef);
				return
            end
			
            if ~ismatrix(arraydef) || size(arraydef,1) < numel(obj.elements)
				[nr,nc,k] = size(arraydef);
				arraydef = reshape(arraydef,nr*nc,k);
				N = nr*nc;
			else
				[N,k] = size(arraydef);
            end
            
            % if N>1 && numel(obj.elements) ~= N, error('Array definition size does not match pvArea elements'); end
            assert(numel(obj.elements)== N,'Array definition size does not match pvArea elements');
           
            if size(unique(arraydef,'rows'),1) < N
				error('At least one element-address appears twice in the array definition')
            end
            if any(arraydef(:)==0)
				warning('reshapeArea:zeroindices','One or more elements will not appear in the new object')
            end
            
            if N == 0
               newobj = pvArea(obj.border);
               return;
            end

            if k == 1
            % just sort existing elements
				newobj = pvArea(obj.border, obj.elements(arraydef));
                return;
            end
            
            % Get a list of valid group indices for the highest hierarchy level of arraydef 
            groupindices = unique(arraydef(:,k));
            validindices = groupindices > 0 & ~isnan(groupindices);
            groupindices = groupindices(validindices);
            if max(groupindices)~= numel(groupindices)
                warning('reshapeArea:missingindex','Object will contain empty elements at level %d',k)
            end

            % Place the corresponding elements of obj inside every new element, by recursively
            % calling reshapeArea with a reduced array definition:
            nwel(max(groupindices),1) = pvArea();
            for j = groupindices'
                inthisgroup = arraydef(:,k) == j;
                nwel(j) = reshapearea(obj.elements(inthisgroup),arraydef(inthisgroup,1:k-1));
            end

            % Create new parent object
            newobj = pvArea(obj.border,nwel);
		end
		
		function arrdf = getarraydefinition(obj,prettyplanes)
		% for a k-depth pvArea object obj with a total of N sub-elements, returns a N·(k-1) array of
        % indices such that obj = reshapearea(flattenpvarea(obj),arraydef). In other words, for 
        % every row of arraydef: a(j,1)...a(j,k-1), if flatobj = flattenpvarea(obj), then:
        %
		%	obj.elements(a(j,k-1)).elements(a(j,k-2))...elements(a(j,1)) = flatobj.elements(j)
        %
        % If prettyplanes = true (default is false), the function will try to return an n·m·k array
        % instead, such that each plane (n,m,j) approximately resembles the physical layout of the 
        % modules. This is meant to work for regular, rectangular arrays (same pitch and azimuth).
		
            if nargin < 2, prettyplanes = false; end
			k = obj.depth;
			S = obj.dims;
			
            if k > 2
				Ne = prod(S)/S(1); % maximum expected number of sub-elements per element
				arrdf = zeros([k-1,Ne,S(1)]);
				used = false(Ne,S(1));
				for j = 1:S(1)
					thisarraydef = getarraydefinition(obj.elements(j),false);  % Nj·(k-2) array, with Nj <= Ne
					used(1:size(thisarraydef,1),j) = true;
					arrdf(1:k-2,used(:,j),j) = shiftdim(thisarraydef,1);	 % store index for lower levels
					arrdf(k-1,used(:,j),j) = j;							 % add index for current level
				end
				arrdf = shiftdim(arrdf(:,used),1);
			else
				arrdf = (1:S)';
            end
            
            if ~prettyplanes 
            % Return table as it is
                return; 
            else
            % Try to find a regular array of indices that resembles the physical object
                
                % Get the centerpoints of all the smallest elements of obj
                flatobj = flattenpvarea(obj);
                Ne = size(arrdf,1);
                idx = zeros(2,Ne);
                for j = 1:Ne
                    [w,h,idx(:,j)] = rectangleproperties(flatobj.elements(j).border);
                end

                % Put them in units of element width and height, to recreate a similar plane pattern
                idx(1,:) = round((idx(1,:)-min(idx(1,:)))/w)+1;
                idx(2,:) = round((max(idx(2,:))-idx(2,:))/h)+1;  % lower elements have higher row index
                idx = bsxfun(@plus,idx,-min(idx,[],2)+1);        % force indices to start at one

                if size(unique(idx','rows'),1) < Ne
                    warning('Could not find a nice pattern to place indices on sheets, you''ll have to do with columns');
                else	
                    % for each hierarchical level, put the indices in a new plane
                    arrtable = arrdf; 
                    arrdf = -ones(max(idx(2,:)),max(idx(1,:)),k-1); % -1 in all unconnected modules.
                    for j = 1:k-1
                        A = sparse(idx(2,:),idx(1,:),arrtable(:,j));% fit in a square array 
                        arrdf(:,:,j) = full(A);
                        %B(~spones(A)) = -1;
                    end
                end
            end
		end

        function P = areaenvelope(A,offset,n)
        % P = AREAENVELOPE(A) - Get a simplified polygonal envelope for the elements of A. i.e.
        %   a simple polygon (single area, no holes) that contains all A.element.border(k)  
        %   polygons for all k in 1:numel(A.elements).
        %
        % P = AREAENVELOPE(...,OFFSET,N) - pass additional parameters to POLYGON.ENVELOPE.
        %
        % See also: ALPHASHAPE, POLYOUTLINE, CONVXHULL, POLYOFFSET
        
            if isempty(A.elements), P = polygon(); return; end
            P = fixorientation([A.elements.border]);
            % P = polyoutline(P,varargin{:});
            
            if nargin < 3 || isempty(offset), offset = 0; end
            if nargin < 4 || isempty(n), n = 6; end

            x = [P.x]';
            y = [P.y]';
            if offset > 0
               p = polygon(n);
               x = x + offset*[0,p.x];
               y = y + offset*[0,p.y];
            end
            k = convhull(x(:),y(:),'simplify',true);
            P = polygon(x(k),y(k));
		end
        
        function plotArea(A,detail,labeling,opacity)
        % plotArea(A,detail,labeling,opacity)
		% Plots a pvArea structure up to depth 'detail', by superimposing semi-transparent patches
		% whose color is determined by their index at any given depth. The result is an image of 
		% the area array in which elements belonging to the same groups (i.e. cells in the same
		% modules, and modules in the same strings) have similar colors.
		
		% if labeling is logical true, then for each of the smallest plotted elements a label with 
		% indices i.j.k... n (where i,j,k... are the parent elements) will be added.
		% if labeling is a string, then all labels will start with this string
		% (this is used for recursive calls)
		% if opacity is specified, the first level or patches will have this given opacity, and
		% further levels will be plotted on top with 0.25·opacity (also used for recursivity).
		
            colorgradient = @hsv;   % constants
            colorbase = 7;
            %layeropacity = 0.2;
			
			if nargin < 2 || isempty(detail), detail = 3; end
			if nargin < 3 || isempty(labeling) || (islogical(labeling)&&~labeling)
				numbered = false; 
            else
				numbered = true;
				if ischar(labeling), parentlabel = labeling; else, parentlabel = char.empty; end
			end
			if nargin < 4, opacity = 1; end
            detail = min(detail,A.depth);
            
            Ne = numel(A.elements);
            if Ne == 1, layeropacity = 0.6; else, layeropacity = 0.2; end
			
            % Colorbase-sorted color list (1,Nc/B, 2·Nc/B... Nc,2,Nc/B+1,2·Nc/B+1,...)
            % provides good contrast among neighboring elements.
            Nc = ceil(Ne/colorbase)*colorbase;
			colorlist = colorgradient(Nc);
            colorlist = colorlist(reshape(reshape(1:Nc,[],colorbase)',[],1),:);
            
            if numbered
                if ~isempty(parentlabel), parentlabel = [parentlabel '.']; end
                labels = arrayfun(@(j) sprintf('%s%d',parentlabel,j),(1:Ne)','unif',0);
            else
                labels = cell(Ne,1);
                labels(:) = {false};
            end
            
            [V,~,F] = poly2vef([A.elements.border],'unif');
            C = colorlist(1:A.dims(1),:);
            patch('Faces',F,'Vertices',V,'FaceVertexCData',C,'EdgeColor','none',...
                    'FaceColor','flat','FaceAlpha',opacity,'EdgeAlpha',opacity,'LineWidth',0.1);
                
%             %plot outline
%             if ~isvoid(A.border)
%                 hold on
%                 plot([A.border.x(:);A.border.x(1)],[A.border.y(:);A.border.y(1)];
%                 hold off
%             end
                
            if detail > 2
                for j = 1:Ne
                    plotArea(A.elements(j),detail-1,labels{j},opacity*layeropacity);
                end
            end

%             for j = 1:Ne
% 				p = patch(A.elements(j).border.x,A.elements(j).border.y,j);
% 				set(p,'FaceColor',colorlist(j,:))
% 				set(p,'FaceAlpha',opacity);
%                 set(p,'LineWidth',0.1);
%                 set(p,'EdgeColor',colorlist(j,:));
%                 if opacity >0, set(p,'EdgeAlpha',opacity); end
% 				axis equal
% 
%                 if detail > 2
% 					plotArea(A.elements(j),detail-1,labels{j},opacity*layeropacity);
%                 else
%                     centers(j,1:2) = centroid(A.elements(j).border)';
%                 end
%             end
            
            if numbered && detail <= 2
                centers = arrayfun(@(e) centroid(e.border)',A.elements,'unif',0);
                centers = cat(1,centers{:});
                centers(:,3) = 0.1;
                text(centers(:,1),centers(:,2),centers(:,3),labels,...
                    'horizontalalignment','center','verticalalignment','middle','fontsize',7);
            end
            axis equal
        end
        
        function plotelements(A,detail,labeling,varargin)
		% Plot the elements of pvArea A at depth 'detail' as patches colored according to their
		% parent elements, mimicking the behaviour of plotArea without the clutter of upper layers
		
		% if labeling is logical true, then for each plotted element a label with 
		% indices i.j.k... n (where i,j,k... are the parent elements) will be added.
		% if labeling is a string, then all labels will start with this string
		% (this is used for recursive calls)
        
        % The two optional arguments in varargin (opacity and parentcolor) are intended for
        % recursive iteration: the color palette for the children elements is set as a variation
		% from color 'parentcolor', with variability controlled by 'opacity'
		
            colorgradient = @hsv;   % constants
            colorbase = 12;
			
            Ne = numel(A.elements);
            if Ne == 0, return; end % Quit if nothing to plot
            
            if Ne == 1, layeropacity = 0.6; else, layeropacity = 0.25; end
                
			if nargin < 2 || isempty(detail), detail = min(3,A.depth); 
            else, detail = min(detail,A.depth); end
            if nargin < 3 || isempty(labeling) || (islogical(labeling)&&~labeling)
				numbered = false; 
            else
				numbered = true;
				if ischar(labeling), parentlabel = labeling; else, parentlabel = char.empty; end
            end
            
            % Start with colorbase-sorted color list (1,Nc/B, 2·Nc/B... Nc,2,Nc/B+1,2·Nc/B+1,...)
            % to provide good contrast among neighboring elements.
            Nc = ceil(Ne/colorbase)*colorbase;
			colorlist = colorgradient(Nc);
            
            colorlist = colorlist(reshape(reshape(1:Nc,[],colorbase)',[],1),:);
            
            if numel(varargin)~=2       % At first function call:
                opacity = 1;          % - use full colormap span
                parentcolor = [1 1 1];  % - start from white (doesn't really matter)
            else
               opacity = varargin{1};
               parentcolor = varargin{2};
               
               % colorlist = colorlist - 0.6*colorlist.*(colorlist./sum(colorlist,2)*parentcolor');
               % colorlist = colorlist + 0.6*(1-mean(colorlist,'all'));
            end
            
            % r = rand(Nc,1)*0.05;
            % colorlist = colorlist.*(1-r) + r;

            % Naïve substractive color mix, should work just fine
            colorlist = 1-bsxfun(@plus,(1-parentcolor)*(1-opacity),(1-colorlist)*opacity);
                        
            % %plot outline
            % hold on
            % plot([A.border.x(:);A.border.x(1)],[A.border.y(:);A.border.y(1)],'r-');
            % hold off
            
            % Create list of labels as parentlabel.elementNo
            if numbered
                if ~isempty(parentlabel), parentlabel = [parentlabel '.']; end
                labels = arrayfun(@(j) sprintf('%s%d',parentlabel,j),(1:Ne)','unif',0);
            elseif detail > 2
            % Pass false to any recursive calls
                labels = cell(Ne,1);
                labels(:) = {false};
            end
            
            % For nested pvAreas, make recursive calls and get out
            if detail > 2
                for j = 1:Ne
					plotelements(A.elements(j),detail-1,labels{j},opacity*layeropacity,colorlist(j,:));
                end
                return;
            end
            
            % Actual plotting: gather info. to make a single call to patch()
            %   - get a list of vertices, and count them
            verts = cell(Ne,1);
            Nv = zeros(Ne,1);
            for j = 1:Ne
                Nv(j) = numel(A.elements(j).border.x);
                verts{j} = [A.elements(j).border.x;A.elements(j).border.y]';
            end
            verts = cat(1,verts{:});
            
            %   - generate connectivity matrix (each row contains a list of vertices for face j)
            faces = nan(Ne,max(Nv));
            for j = 1:Ne
                faces(j,1:Nv(j)) = (1:Nv(j)) + sum(Nv(1:j-1)); 
            end
  
            %   - reduce colorlist to Ne faces
            colorlist = colorlist(1:Ne,:);

            patch('Faces',faces,'Vertices',verts,'FaceVertexCData',colorlist,...
                            'FaceColor','flat','LineWidth',0.1,'EdgeColor',[0.8,0.8,0.8]);

            % Get label positions, and plot using a single call to text()
            if numbered
                centers = zeros(numel(A.elements),3);
                centers(:,3) = 0.1;
                for j = 1:Ne
                    centers(j,1:2) = centroid(A.elements(j).border)';
%                     centers(j,1) = 0.5*(min(A.elements(j).border.x)+max(A.elements(j).border.x));
%                     centers(j,2) = 0.5*(min(A.elements(j).border.y)+max(A.elements(j).border.y));
                end
                text(centers(:,1),centers(:,2),centers(:,3),labels,...
                    'horizontalalignment','center','verticalalignment','middle','fontsize',7);
            end
            %axis equal0
		end
        
		function newobj = reducedetail(obj,maxdetail)
		% eliminates any sub-elements beyond depth 'maxdetail'
		% if maxdetail is ommited, the last level of detail is taken away

			if nargin < 2, maxdetail = obj.depth-1; end
			
            newobj = copy(obj); % start with a copy	
            
			% rule out trivial case in which there's nothing to be done
            if maxdetail <= 0 || obj.depth <= maxdetail, return; end
				
			if maxdetail == 1
            % desired level reached, delete all elements below
                newobj.elements = pvArea.empty;
            else
            % apply recursively to obj elements, with maxdetail-1
                newobj.elements = arrayfun(@(e) reducedetail(e,maxdetail-1),obj.elements);
			end
		end

        function S = struct(A,varargin)
        % S = STRUCT(A) - Recursively convert a pvArea object to a (nested) structure
        % s = STRUCT(A,'pack[16]',[SCALE]) - pack/pack16 every border polygon
        
            % MATLAB seems too stupid to tell appart struct(A,..) and struct(..,'field',A)
            if ~isa(A,'pvArea')
                S = builtin('struct',A,varargin{:});
                return;
            end
            
            [opt,rem] = getflagoptions(varargin,{'pack','pack16'});
            if isempty(rem), opt.scale = [];
            else
               assert(numel(rem) == 1 && isempty(rem{1}) || isscalar(rem{1}),'Unexpected arguments')
               opt.scale = rem{1};
            end
            S = recursivestruct(A,opt);

            function S = recursivestruct(A,opt)
                if isempty(A)
                    pp = {'border','elements','area','depth','dims'};
                    S = cell2struct(cell(numel(pp),0),pp);
                    return;
                end
                for j = numel(A):-1:1
                    if opt.pack16
                        S(j).border = pack16(A(j).border,opt.scale);
                        S(j).area = polygon.packedarea(S(j).border);
                    elseif opt.pack
                        S(j).border = pack(A(j).border,opt.scale);
                        S(j).area = polygon.packedarea(S(j).border);
                    else
                        S(j).border = A(j).border;
                        S(j).area = A(j).area;
                    end
                    S(j).depth = A(j).depth;
                    S(j).dims = A(j).dims;
                    if ~isempty(A(j).elements)
                        S(j).elements = recursivestruct(A(j).elements,opt);
                    else
                        S(j).elements = struct.empty;
                    end
                end
                S = reshape(S,size(A));
            end
        end

        function c = firstleaf(A)
        % C = FIRSTLEAF(A) - returns the first smallest sub-element of a pvArea object A.
        % That is: C = A.element(1).element(1)...element(1), such that C.elements is empty.

            c = A;
            switch numel(A)
            case 0, return;
            case 1
                while ~isempty(c.elements), c = c.elements(1); end
                c = recenter(copy(c));
            otherwise
                error('FIRSTLEAF only works for scalar objects');
            end
            
            function c = recenter(c)
                r = centroid(c.border);
                c.border = polytranslate(c.border,-r);
            end
		end

    end
    methods (Static = true)
        function obj = loadobj(S)      
            if isempty(S), obj = pvArea.empty(); return; end
            if isa(S,'pvArea') && ~isempty(S.area), obj = S; return; end  % already loaded
            
            if isstruct(S.border) && isfield(S.border,'scale') && isinteger([S.border.x])
            % Unpack border polygons, if required
                S.border = polygon.unpack(S.border);
            end
            obj = pvArea(S.border,arrayfun(@pvArea.loadobj,S.elements));
        end
       function [ShF,varargout] = shadingfactors(obj,P,varargin)
        % [SHF,h] = SHADINGFACTORS(OBJ,P,...)
        % Calculates the percentage coverage of a polygon(s) P over each element of obj.
        % INPUT:
        %   OBJ - PVAREA object
        %   P - polygon, packed polygon, or cell-array of such
        %   ..'tol',F - fraction of smallest subelement of OBJ that can be disregarded for
        %       simplification, assuming that the whole area is divided into sub-children.
        %       i.e. ShF(:) = 0 as long as the area of intersection of P and obj.border is less 
        %       than obj.border.area/prod(obj.dims)·tol, and similarly for ShF(:) = 1.
        %       Default is SimOptions.cellshadingthreshold
        %   ..'-pack16' - Use pack16(p) instead of pack(p) when packing polygons for clipping
        %   ..'depth',D - Set max. depth level of recursive calculation. 
        %       Equivalent to PVAREA.SHADINGFACTORS(REDUCEDETAIL(OBJ,D),..)
        %       NOTE: the threshold for simplification is still obj.border.area/prod(obj.dims)·tol
        %       and not obj.border.area/prod(obj.dims(1:D-1))·tol
        %   ..'scale',S - Use scale S when packing polygons for clipping
        %   ..'-plot' - plot shaded-element patches on current figure.
        %   ..'-color',c,'edgecolor',e,'alpha',a - parameters for patch plotting
        %
        % OUTPUT:
		%   ShF - an array of size obj.dims (or [numel(P),obj.dims]) with area coverage fractions.
		%     If shadepolygon actually represents a shadow, ShF(j) = 0, means no shade over ShF(j)
		%     If shadepolygon is a light polygon, ShF(j) = 0, means no light, ShF(j) = 1 no shade
        %   h - if '-plot' flag is used, returns vector of plotted patch handles.
                          
            [opt,varargin] = getflagoptions(varargin,{'-plot','-pack16'});
            opt.tol = getSimOption('cellshadingthreshold');
            opt.scale = [];
            opt.depth = obj.depth;
            opt.color = 'k';
            opt.edgecolor = 'k';
            opt.alpha = 1/2;
            opt.waitbar = ~opt.plot && numel(P) > 100;
            opt = getpairedoptions(varargin,opt,'restchk');
            
            opt.depth = min(opt.depth,obj.depth);

            % pack pvArea object
            if isa(obj,'pvArea')
                if opt.pack16, flag = 'pack16'; else, flag = 'pack'; end
                obj = struct(obj,flag,opt.scale);
            else
                assert(isstruct(obj) && all(isfield(obj,{'dims','area','border','elements'})),...
                    'OBJ does not look like a [packed] PVAREA object');
            end
            
            if opt.plot
                h = matlab.graphics.primitive.Patch.empty;
                varargout{1} = h;
            end
            
            if isa(P,'polygon') || isstruct(P) && all(isfield(P,{'x','y','scale'}))
            % single polygon calculation
                ShF = zeros([obj.dims,1]);
                if isempty(P), return; end
                ShF = recursivefactors(obj,P,opt);
            else
            % cell-array of polygons
                assert(iscell(P),'P must be a [packed] polygon or cell-array of poylgons');
                ShF = zeros([numel(P),obj.dims(1:opt.depth-1),1]);
                if isempty(P), return; end
                
                if opt.waitbar
                    wb = optwaitbar(0,'Calculating shading factors...');
                end
                
                for k = 1:numel(P)
                    if opt.plot, delete(h); h(:) = []; end
                    fk = recursivefactors(obj,P{k},opt);
                    ShF(k,:) = fk(:);
                    if opt.plot, drawnow(); end
                    if opt.waitbar
                        wb.update(k/numel(P),'Calculating shading factors...','-addtime'); 
                    end
                end
            end
            if opt.plot, varargout{1} = h; end

            function ShF = recursivefactors(G,P,opt)
                
                if ~isstruct(P)  % pack polygon, if not already
                    if opt.pack16, P = pack16(P,opt.scale); else, P = pack(P,opt.scale); end
                end
                
                thresh = opt.tol/prod(G.dims); % threshold is relative to full-depth!
                ShF = zeros([G.dims(1:opt.depth-1),1]);
                [f,P] = shf(P,G.border,G.area);
                if isempty(P) || f < thresh, return; end
                if f > 1-thresh, f = 1; end
                
                if f < 1 && ~isempty(G.elements) && opt.depth > 1
                    
                    if opt.plot, last_plotted = numel(h); end

                    opt.depth = opt.depth - 1;
                    for j = 1:G.dims(1)
                       ff = recursivefactors(G.elements(j),P,opt);
                       ShF(j,:) = ff(:);
                    end
                    
                    if all(ShF == 1,'all') && opt.plot
                        f = 1;
                        delete(h(last_plotted+1:end));
                        h(last_plotted+1:end) = [];
                    end
                end
                
                if f == 1 || isempty(G.elements) || opt.depth == 1
                    ShF(:) = f;
                    if opt.plot
                        upck = @(x) double(x)/double(G.border.scale);
                        h(end+1) = patch(upck(G.border.x),upck(G.border.y),opt.color,...
                                    'Facealpha',f*opt.alpha,'EdgeColor',opt.edgecolor); 
                    end
                else
%                     opt.depth = opt.depth - 1;
%                     for j = 1:G.dims(1)
%                        f = recursivefactors(G.elements(j),P,opt);
%                        ShF(j,:) = f(:);
%                     end
                end
            end
            
            function [f,P] = shf(P,G,a)
            % Fraction WITHOUT light-polygon LP on pvArea: f = [1-area(G & P)/area(G)]·MAX_F

                P = polyclip(G,P,1); % G.border & LP
                f = polygon.packedarea(P)/a;
            end
        end
    end
end
