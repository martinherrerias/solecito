classdef polygon
% POLYGON class - TODO: should probably be replaced by POLYSHAPE
%
% Represents flat, closed, polygons as structure-arrays with fields x, y, and hole. A simple 
% poligon with n vertices has size(p.x) = size(p.y) = [1,n], i.e. the first vertex is assumed to
% also be the last. Complex polygons (e.g. polygons with holes) are represented as arrays of simple
% polygons, some of which can be marked as holes.
%
% NOTE: for some functions (e.g. clipping operations) vertex order can affect the way holes are 
% interpreted (see fill-types in POLYCLIP). The polygon constructor ensures that vertex-order
% works with most fill-types, and the default fill-type is Even-Odd, so that 'true' holes
% (i.e. polygons that fit in a positive polygon of the same array) remain holes even after 
% vertex-order manipulation. Use FIXORIENTATION and check for SELFINTERSECT if things behave
% strangely.
%
% PROPERTIES
%   x: horizontal vector of x-coordinates
%   y: horizontal vector of y-coordinates
%   hole: boolean value
%   area (Dependent) - use POLYGONAREA for complex polygons
%
% METHODS

    properties
        x       % horizontal vector of x-coordinates
        y       % horizontal vector of y-coordinates
		hole    % boolean scalar
	end
	properties (Dependent = true)
        signedarea
        perimeter
    end
    properties (Constant)
       SCALE = int64(2)^32    % ~2e-10 precision, room for vertex coordinates within ± 1e9
       SCALE16 = int16(2)^10  % ~1e-3 precision, room for vertex coordinates within ± 32 
    end
    methods
        function obj = polygon(varargin)
        % POLYGON class constructor:
		% POLYGON()	- Empty polygon (1x1 polygon array with empty fields)
        %       Use POLYGON.EMPTY() for a 0x0 array of polygons.
        % POLYGON(n) - Regular polygon of n sides, with unit radius.
		% POLYGON(S) - Polygon from [array of] structure(s). S(j).x, S(j).y are vectors with
		%       the vertices of a poligon j. Optionally S(j).hole defines whether S(j) is
		%       a cut-out region.
        % POLYGON(xx,yy) - Creates a simple POLYGON object whose vertices are defined by pairs
        %       xx(j), yy(j) for j = 1:numel(xx).
        % POLYGON(xy) - Polygon from 2xN array of vertices xy
		% POLYGON(w,h) - creates a rectangle of w·h, centered at origin
		% POLYGON([x1,x2],[y1,y2]) - creates a rectangle within the given boundaries
		% POLYGON(TH,R,'pol') - defines a poligon in polar coordinates. Angles in TH must be
        %       provided in degrees, R can be a scalar, and will be replicated for all TH.
		% POLYGON(..,..,hole) - hole (boolean) defines if the polygon is a cut-out region.
        %
        % See also: POLYGON3D, POLYGON
		
            narginchk(0,4);
            warning_reseter = naptime();  %#ok<NASGU>
            warning('off','polygon:areasign');
            
            ishole = false;
            if ~isempty(varargin) && islogical(varargin{end}) || ...
                    (nargin > 2 && isnumeric(varargin{end}))
            % if last argument is 'hole', include in polygon definition, then dump
                ishole = varargin{end} > 0;
                varargin(end) = [];
            end
                
            switch numel(varargin)
            case 0
            % polygon() - empty polygon constructor
                obj.x = [];
                obj.y = [];
                obj.hole = ishole;
            case 1
                if isempty(varargin{1})
                % from an empty array of structures, return an empty array of polygons
                    obj = polygon.empty;
                    return;
                elseif isstruct(varargin{1}) || isa(varargin{1},'polygon')
                % polygon(S) or polygon(P3D)
                    P = varargin{1};
                    if isstruct(P) && isfield(P,'scale') && isinteger([P.x])
                       warning('Assuming packed polygon. Use UNPACK or remove P.scale field');
                       P = polygon.unpack(P);
                    end
                    
                    for j = numel(P):-1:1
                        obj(j).x = P(j).x(:)';
                        obj(j).y = P(j).y(:)';
                    end

                    if ~isstruct(P) || isfield(P,'hole')
                        [obj.hole] = deal(P.hole); 
                    else
                        [obj.hole] = deal(ishole); 
                    end
                    obj = fixorientation(obj);
                    return
                elseif isscalar(varargin{1})
                % polygon(n)
                    n = varargin{1};
                    assert(mod(n,1)==0 && n > 2,'polygon(n) requires an int > 2 number of sides');
                    th = (1/2:1:n)*360/n - 180;
                    obj = polygon(th,ones(size(th)),'pol');
                    obj.hole = ishole;
                    obj = fixorientation(obj);
                    return
                else
                % polygon(xy)
                    xy = varargin{1};
                    if size(xy,1)~=2, xy = xy'; end
                    assert(isnumeric(xy) && size(xy,1)==2,...
                        'polygon:narg1','Expecting structure(s), polygon3d or XY array');
                    obj.x = xy(1,:);
                    obj.y = xy(2,:);
                    obj.hole = ishole;
                    obj = fixorientation(obj);
                    return
                end
                otherwise
                if numel(varargin) > 2
                    assert(ischar(varargin{3}) && strcmpi(varargin{3},'pol'),...
                                        'Expecting ''pol'' as third argument');
                    % transform from polar coordinates, if required 
                    [varargin{1},varargin{2}] = pol2cart(varargin{1}*pi/180,varargin{2});
                end

                xx = varargin{1};
                yy = varargin{2};
                if numel(xx) ~= numel(yy)
                    error('polygon:narg2','Not sure how to read the two arguments');
                end

                %if numel(varargin) == 2
                    switch numel(xx)
                    case 1 
                        % create x·y rectangle centered at 0,0 (clockwise vertices)
                        xx = [-0.5*xx 0.5*xx 0.5*xx -0.5*xx];
                        yy = [-0.5*yy -0.5*yy 0.5*yy 0.5*yy];
                    case 2
                        % create rectangle from limits x(1:2) y(1:2)
                        xx = [xx(1) xx(2) xx(2) xx(1)];
                        yy = [yy(1) yy(1) yy(2) yy(2)];
                    % otherwise leave xx,yy vertices as they are
                    end
                %end

                % store whatever polygon x,y represents as element 1
                obj.x = xx(:)';
                obj.y = yy(:)';
                obj.hole = ishole;
                obj = fixorientation(obj); % set vertex orientation according to obj.hole
            end
        end
		
        % Orientation
        function isccw = orientation(p)
        % isccw = ORIENTATION(p) determine if the vertices of polygon(s) p are oriented counter-
        % clock-wise when the normal faces up, i.e. around the normal following the right-hand-rule
        
            [p.hole] = deal([]);    % remove hole-flags to avoid warnings
            isccw = reshape([p.signedarea] >= 0,size(p));
        end
        function p = reverse(p,holesalso)
        % flip vertex order, so that all(orientation(reverse(p)) = !orientation(p))
        % if holesalso, flip also the p.hole flag of all elements of p
            if nargin < 2, holesalso = true; end
            
            for j = 1:numel(p)
               p(j).x = fliplr(p(j).x); 
               p(j).y = fliplr(p(j).y);
               if holesalso, p(j).hole = ~(p(j).hole); end
           end
        end
        function p = fixorientation(p,fixholes)
        % Q = FIXORIENTATION(P,FIXHOLES) - Check that all 'hole' flags in P match vertex winding 
        %   convention CCW+. When this is not the case:
        %       + if FIXHOLES == false (default), vertices will be re-ordered to match p.hole
        %       + if FIXHOLES == true, hole-flags will be set to match vertex orientation
           
            if isvoid(p), p = polygon.empty(); end
            if nargin < 2, fixholes = false; end
                            
            wrong = ~xor(orientation(p(:)),[p.hole]');
            if ~any(wrong), return; end
            if fixholes
                for j = find(wrong)', p(j).hole = ~p(j).hole; end
            else
                for j = find(wrong)'
                   p(j).x = fliplr(p(j).x); 
                   p(j).y = fliplr(p(j).y);
                end
            end
        end
        
        % Polygon propperty query
		function yn = isvoid(P)
        % yn = ISVOID(P) - return true for an empty polygon array P, or when all elements of 
        % polygon array P have less than three vertices. NOTE that ~isvoid(P) does not guarantee
        % that area(P) > 0 (vertices can be redundant, or holes larger than positive areas).
        %
        % See also: POLYGONAREA, POLYGON
            yn = true;
			if isempty(P), return; end
            for j = 1:numel(P)
                if numel(P(j).x) > 2, yn = false; break; end
            end
        end    
        function A = get.signedarea(p)
		% Returns the area of individual polygons using the shoelace-formula
        
            if isempty(p), A = 0; return; end       
            if numel(p.x) < 3; A = 0; return; end % Degenerate polygon
            
            cx = (min(p.x)+max(p.x))/2;
            cy = (min(p.y)+max(p.y))/2;
            
            xi = p.x - cx;   % reduce |x| and |y| to increase numerical precision
            yi = p.y - cy;
            xj = circshift(xi,-1);  % faster than x([2:end,1])
            yj = circshift(yi,-1); 
                        
            % A = sum(xi·yj - yj·xi)/2 = sum[(xi-xj)·yj - (yi-yj)·xj]/2, better precision
            
            T1 = (xi-xj).*yj;
            T2 = (yi-yj).*xj;
            A = sum(T1 - T2)/2;
            
            tol = max(eps([T1,T2]));
            if abs(A) <= tol, A = 0; return; end            
            
            if ~isempty(p.hole) && abs(A) > eps && sign(A) ~= 1-2*p.hole
               warning('polygon:areasign','Inconsistent winding & polygon.hole'); 
            end
        end
        function A = area(P), A = sum([P.signedarea]); end
        function S = get.perimeter(p)
            S = sum(sqrt((p.x - circshift(p.x,-1)).^2 + (p.y - circshift(p.y,-1)).^2));
        end
        function S = length(P), S = sum([P.perimeter]); end
			
        function c = centroid(p)
        % C = CENTROID(P) return the centroid [2,1] vector of (complex) polygon P.
        %   xc = sum((xi + xj)·ti)/6 where ti = xi·yj - yj·xi
            a = zeros(numel(p),1); xa = a; ya = a;
            for j = 1:numel(p)
                x0 = mean(p(j).x);  % reduce |x| and |y| to increase numerical precision
                y0 = mean(p(j).y);
                xi = p(j).x - x0;
                yi = p(j).y - y0;
                xj = circshift(xi,-1);  % faster than x([2:end,1])
                yj = circshift(yi,-1); 
                t = (xi-xj).*yj - (yi-yj).*xj; % == xi·yj - xj·yi
                a(j) = sum(t)/2;
                xa(j) = sum((xi + xj).*t)/6 + x0*a(j); % reapply offset x0, y0
                ya(j) = sum((yi + yj).*t)/6 + y0*a(j);
            end
            c = [sum(xa);sum(ya)]/sum(a);
        end
        function I = secondmoment(p,q)
        % I = SECONDMOMENT(P,Q) - Return second area moment of polygon P around point Q in matrix
        %   form: [Ixx,Ixy;Iyx,Iyy]. If Q is ommited, default is centroid(P)
            if nargin < 2, q = centroid(p); end
            xx = 0; yy = 0; xy = 0;
            for j = 1:numel(p)
                xi = p(j).x - q(1); 
                yi = p(j).y - q(2);
                xj = xi([2:end,1]);
                yj = yi([2:end,1]); 
                t = xi.*yj - xj.*yi;
                xx = xx + sum((xi.^2 + xi.*xj + xj.^2).*t)/12;
                yy = yy + sum((yi.^2 + yi.*yj + yj.^2).*t)/12;
                xy = xy + sum((xi.*yj + 2*xi.*yi + 2*xj.*yj + xj.*yi).*t)/24;
            end
            I = sign(xx)*[xx,xy;xy,yy];
        end
        function [c,r] = maxinscribedcircle(P,q)
        % [c,r] = MAXINSCRIBEDCIRCLE(P) - return the APPROXIMATE center and radius of the maximum
        %   inscribed circle in (complex) polygon P.
        % [c,r] = MAXINSCRIBEDCIRCLE(P,Q) - include additional point(s) Q to the maximization
        %   problem, e.g. to find the best place for a new label considering existing ones.
        %
        % INPUT:    P - POLYGON object, can be complex (i.e. array of polygons)
        %           Q - N·2 array of additional points
        % OUTPUT:   c - 1·2 array, center of maximum inscribed circle
        %           r - radius of maximum inscribed circle.
        %
        % PROVISIONAL: currently uses a quick & dirty approx. should be replaced by a robust,
        %   more precise algorithm.
        
            DS = sum([P.perimeter])/100; 
            
            if nargin < 2, q = zeros(0,2); end
            if size(q,2)~=2,q = q'; end
            assert(size(q,2) == 2,'Additional points must be provided in Nx2 array');

            if isvoid(P), c = zeros(0,2); r = NaN; return; end
            
            p = refinepoly(P,DS);        % add intermediate points to approximate edge constraints
            V = poly2vef(p);             % list all polygon vertices
            V = cat(1,q,V);              % add any previously chosen points as vertices
            u = voronoin(V);             % find Voronoi-cell-vertices inside polygon
            u = u(insidepolygon(P,u(:,1),u(:,2)),:);
            D = min(pdist2(V,u),[],1);   % distance from each Voronoi point to closest vertex/edge
            [r,j] = max(D);
            c = u(j,:);
        end
		function [w,h,c,d] = rectangleproperties(P)
        % [w,h,c,d] = RECTANGLEPROPERTIES(P)
		% Returns the maximum span in x and y of polygon P, along with the centerpoint and length  
        % of the diagonal for the rectangle that fits around the polygon.

			w = (max(P.x)-min(P.x));
			h = (max(P.y)-min(P.y));
			c = zeros(2,1);
			c(1) = (max(P.x)+min(P.x))/2;
			c(2) = (max(P.y)+min(P.y))/2;
            
            d = hypot(w,h);
        end
        function [xy,idx] = selfintersect(p)
        % [XY,IDX] = SELFINTERSECT(P) - find any (self) intersections of a (set of) polygon(s) P.
        %
        %   P: polygon, or array of polygons.
        %  XY: n·2 array of coordinates for intersection points.
        % IDX: n·4 array of indices locating the intersecting segments, where each row j
        %       IDX(j,:) = [A,B,C,D] implies that edge A in P(B) crosses edge C in P(D).
        %       NOTE that edge k in a polygon q joins vertices [q,q+1], or [n,1] when k = n.       
        %
        % Equations (*) borrowed from: http://algs4.cs.princeton.edu/91primitives/
        % FUTURE: This is a complete-search algorithm, meant for a small number of vertices, 
        % it should eventually be replaced by an implementation of the Bentley-Ottman algorithm.
        %
        % See also: POLYGON
            
            [V,E] = poly2vef(p); % get vertex and edges
            
            vx = find(accumarray(E(:),1) > 2); % vertex-crossings
            if ~isempty(vx), xy = V(vx,:); else, xy = zeros(0,2); end
            
            Ne = size(E,1);
            idx = zeros(nchoosek(Ne,2),2);  % idx = nchoosek(1:Ne,2), just way faster
            k = 0;
            for j = 1:Ne
                idx(k+(1:Ne-j),1) = j;
                idx(k+(1:Ne-j),2) = j+1:Ne;
                k = k+Ne-j;
            end
            
            C = [E(idx(:,1),:),E(idx(:,2),:)];
            C = (C(:,1) ~= C(:,3)) & (C(:,1) ~= C(:,4)) & (C(:,2) ~= C(:,3)) & (C(:,2) ~= C(:,4));
            idx = idx(C,:); % combinations of edges without common vertices
            
            % get lists of coordinates for edge-pairs a-b, c-d
            ax = V(E(idx(:,1),1),1);
            ay = V(E(idx(:,1),1),2);
            bx = V(E(idx(:,1),2),1);
            by = V(E(idx(:,1),2),2);
            cx = V(E(idx(:,2),1),1);
            cy = V(E(idx(:,2),1),2);
            dx = V(E(idx(:,2),2),1);
            dy = V(E(idx(:,2),2),2);
            
            % calculate edge-pair intersection points...
            d = ((bx-ax).*(dy-cy)-(by-ay).*(dx-cx));
            r = ((ay-cy).*(dx-cx)-(ax-cx).*(dy-cy))./d;
            s = ((ay-cy).*(bx-ax)-(ax-cx).*(by-ay))./d;
            
            % the edge-pair actually intersects only if 0 <= r <= 1 and 0 <= s <= 1
            f = (r >= 0) & (r <= 1) & (s >= 0) & (s <= 1);
            
            if nnz(f) > 0
                xy = cat(1,xy,[ax(f) + r(f).*(bx(f)-ax(f)), ay(f) + r(f).*(by(f)-ay(f))]);
            end
            
            if ~isempty(xy) && nargout > 1
            % Generate a local edge index: ie(j,:) = [k,m] -> row j in E is edge k in polygon m
                for j = numel(p):-1:1
                    nv = numel(p(j).x);
                    ie{j} = [j(ones(nv,1)),(1:nv)'];
                end
                ie = cat(1,ie{:});
                
                % Return a list idx(j,:) = [a,b,c,d] -> edge a in P(b) crosses edge c in P(d)
                idx = idx(f,:);
                idx = [ie(idx(:,1),:),ie(idx(:,2),:)];
            end
        end
		function in = insidepolygon(P,x,y,FT)
		% in = INSIDEPOLYGON(P,xq,yq) - Generalization of INPOLYGON(), to work directly with 
        %   (compound) polygon objects.
        %
        % in = INSIDEPOLYGON(P,xq,yq,FT) - Specify fill-type (for complex polygons), select
        %       0:3/{'EO','NZ',POS','NEG'} for Even-Odd (default), Non-Zero, Positive,
        %       or Negative. See: http://www.angusj.com/delphi/clipper/documentation/...
        %                         Docs/Units/ClipperLib/Types/PolyFillType.htm
        %
        % See also: INPOLYGON
        
            in = false(size(x));
            if isempty(P), return; end
            
            narginchk(3,4);
            assert(all(size(x) == size(y)),'X,Y must be equal-sized sets of coordinates'); 
            
            if nargin < 4 || isempty(FT), FT = 0; end
            if ischar(FT), FT = find(strcmpi(FT,{'eo','nz','pos','neg'})) - 1; end
            assert(isscalar(FT) && any(FT == 0:3),'polyclip:FT','Unexpected Fill-Type');
    
            % get winding number
            w = zeros(numel(x),1,'int8');
            for j = 1:numel(P)
				if P(j).hole, w = w - int8(inpolygon(x(:),y(:),P(j).x,P(j).y));
                else, w = w + int8(inpolygon(x(:),y(:),P(j).x,P(j).y));
				end
            end

            switch FT
                case 0, in(:) = mod(w,2) ~= 0;
                case 1, in(:) = w ~= 0;
                case 2, in(:) = w > 0;
                case 3, in(:) = w < 0;
            end
        end
        function Q = convhull(P)
        % Q = CONVHULL(P) - Return the convex-hull of (a set of) polygon(s) P
        % See also: CONVHULL
        
            xy = [[P.x];[P.y]]';
            K = convhull(xy);    % convex hull point indices
            Q = polygon(xy(K,:));
            
        end
        
        % Conversion
		function S = poly2struct(P)
		% Convert array of polygons to array of structures
            if isempty(P)
                S = struct('x',{},'y',{},'hole',{});
            else
                for j = numel(P):-1:1
                    S(j).x = P(j).x;
                    S(j).y = P(j).y;
                    S(j).hole = P(j).hole;
                end
            end
        end
        function [V,E,F] = poly2vef(p,unif)
        % [V,E,F] = POLY2VEF(P,UNIF) - Generates a list of vertices, edges, and faces of a (set of) 
        % polygon(s) P.
        %
        %   V: Nv·2 array of vertices, where V(j,:) = [xj,yj]. When numel(p) > 1, the vertices
        %      of all polygons p(k) will be concatenated in order.
        %   E: Ne·2 array of edges (vertex-index-pairs), where E(k,:) = [i,j] implies that there is 
        %      a vertex from V(i,:) to V(j,:).
        %   F: When UNIF = false (default), F is a cell-array of vertex-index-tuples. Each element 
        %      k of F represents a face formed by vertices F{k} = [a,b,..]
        %      When UNIF = true, F is packed into a uniform Nf·t array of (vertex-index-tuples), 
        %      where t is the maximum number of vertices for all faces, and face k is given by 
        %      F(k,~isnan(F(k,:)) i.e. for faces with less than t vertices, non-used positions are
        %      filled with NaN.
        
            if nargin < 2, unif = false; end
        
            Np = numel(p);
            V = cell(Np,1);
            Nv = zeros(Np,1);
            for j = 1:Np
                V{j} = [p(j).x;p(j).y]';
                Nv(j) = numel(p(j).x);
            end
            
            if isempty(Nv) || all(Nv == 0)
                V = zeros(0,2); E = zeros(0,2);
                if unif, F = []; else, F = {}; end
                return;
            end
            
            V = cat(1,V{:});
            [V,~,ia] = unique(V,'rows','stable'); % Remove duplicate vertices
            
            if nargout > 1
                for j = Np:-1:1, E{j} = ia([1:Nv(j);[2:Nv(j),1]]'+sum(Nv(1:j-1))); end
                E = cat(1,E{:});
            end
            
            if nargout > 2
                if ~unif
                    for j = Np:-1:1, F{j} = ia((1:Nv(j)) + sum(Nv(1:j-1)))'; end
                else
                    F = nan(Np,max(Nv));
                    for j = 1:Np
                        F(j,1:Nv(j)) = ia((1:Nv(j)) + sum(Nv(1:j-1))); 
                    end
                end
            end
        end
        function [V,E,F,Issues] = triangulate(P)
        % TRIANGULATE converts polygon p(s) into a set of triangles, using constrainded
        % Dalaunay-triangulation.
        %   p: polygon or array of non-intersecting polygons (for polygons with holes). 
        %      CAUTION: the function isInterior used to determine what is and what's not a hole
        %      doesn't care about p.hole flags, but follows the even-odd fill-rule. Any output of
        %      polygon operations should do the trick, so try p = mergepolygons(p,polygon()) when
        %      in doubt.    

            [V,E] = poly2vef(P);
            F = delaunayTriangulation(V,E);

            if size(F.Points,1)~= size(V,1) %  Intersecting edge constraints
               Issues.V = F.Points(size(V,1)+1:end,:);
               Issues.E = F.Constraints(size(E,1)+1:end,:);
               V = F.Points;
            else
               Issues = struct('V',[],'E',[]);
            end
            
            F = F(isInterior(F),:);
        end
        
       % Function calls to polygon clipping library
        function S = pack(p,scale,hiRange)
        % S = PACK(P) convert polygon P coordinates to int62 values, compatible with CLIPPER
        % PACK is called by default inside POLYCLIP. However, for repetitive operations over the 
        % same set of polygons, using POLYCLIP over ready-packed polygons can be much faster.
        % See also: POLYCLIP, CLIPPER, POLYGON.UNPACK
            if nargin < 2 || isempty(scale), scale = polygon.SCALE; end
            if isvoid(p), S = struct('x',int64(0),'y',int64(0),'scale',scale); return; end
            if nargin < 3, hiRange = int64(2)^62-1; end
            s = double(scale);
            for j = numel(p):-1:1
                S(j).x = min(hiRange,max(-hiRange,int64(p(j).x*s)));
                S(j).y = min(hiRange,max(-hiRange,int64(p(j).y*s)));
                S(j).scale = scale;
            end
        end
        function S = pack16(p,scale,hiRange)
        % S = PACK(P) convert polygon P coordinates to int16 values, 
        % See also: POLYGON.PACK, CLIPPER, POLYGON.UNPACK
            if nargin < 2 || isempty(scale), scale = polygon.SCALE16; end
            if isvoid(p), S = struct('x',int16(0),'y',int16(0),'scale',scale); return; end
            if nargin < 3, hiRange = 2^15-1; end
            s = double(scale);
            for j = numel(p):-1:1
                S(j).x = min(hiRange,max(-hiRange,int16(p(j).x*s)));
                S(j).y = min(hiRange,max(-hiRange,int16(p(j).y*s)));
                S(j).scale = scale;
            end
        end
        function M = mergepolygons(P,Q,varargin)
		% M = MERGEPOLYGONS(P,Q) - returns the union of (complex) polygons P and Q
        % M = MERGEPOLYGONS(P) - returns the union of a set of simple (no-holes!) polygons P
        % M = MERGEPOLYGONS(.. pFT, qFT, scale) passes additional parameters to POLYCLIP.
        %
        % P and Q are expected to be non-intersecting polygon objects, or array of polygon objects
        % representing a complex polygon. Make sure to FIXORIENTATION if P, Q are not the output 
        % of another POLYCLIP function.
        %
        % See also: POLYCLIP, FIXORIENTATION, SUBSTRACTPOLYGONS, INTERSECTPOLYGONS, XORPOLYGONS
        
			if nargin > 1 && isa(Q,'polygon')
            % M = MERGEPOLYGONS(P,Q)
                M = polyclip(P,Q,3,varargin{:});
            else
            % M = MERGEPOLYGONS(P)
                if nargin > 1, varargin = [{Q} varargin]; end
                holes = [P.hole];
                if any(holes) && ~isscalar(P)
                    warning('mergepolygons:arrayholes',['Array-merging a set of polygons '...
                        'that contains holes. Use one-by-one merge for predictable results.']);
                end
                M = polyclip(P,polygon(),3, varargin{:});
			end
        end
		function PQ = substractpolygons(p,q,varargin)
        % M = SUBSTRACTPOLYGONS(P,Q) - returns the difference of (complex) polygons P and Q
        % M = SUBSTRACTPOLYGONS(.. pFT, qFT, scale) passes additional parameters to POLYCLIP.
        %
        % P and Q are expected to be non-intersecting polygon objects, or array of polygon objects
        % representing a complex polygon. Make sure to FIXORIENTATION if P, Q are not the output 
        % of another POLYCLIP function.
        %
        % See also: POLYCLIP, FIXORIENTATION, MERGEPOLYGONS, INTERSECTPOLYGONS, XORPOLYGONS
            PQ = polyclip(p,q,0,varargin{:});
        end
		function PQ = intersectpolygons(p,q,varargin)
        % M = INTERSECTPOLYGONS(P,Q) - returns the intersection of (complex) polygons P and Q
        % M = INTERSECTPOLYGONS(.. pFT, qFT, scale) passes additional parameters to POLYCLIP.
        %
        % P and Q are expected to be non-intersecting polygon objects, or array of polygon objects
        % representing a complex polygon. Make sure to FIXORIENTATION if P, Q are not the output 
        % of another POLYCLIP function.
        %
        % See also: POLYCLIP, FIXORIENTATION, MERGEPOLYGONS, SUBSTRACTPOLYGONS, XORPOLYGONS
            PQ = polyclip(p,q,1,varargin{:});
        end
		function PQ = xorpolygons(p,q,varargin)
		% M = XORPOLYGONS(P,Q) - returns the XORion of (complex) polygons P and Q
        % M = XORPOLYGONS(.. pFT, qFT, scale) passes additional parameters to POLYCLIP.
        %
        % P and Q are expected to be non-intersecting polygon objects, or array of polygon objects
        % representing a complex polygon. Make sure to FIXORIENTATION if P, Q are not the output 
        % of another POLYCLIP function.
        %
        % See also: POLYCLIP, FIXORIENTATION, MERGEPOLYGONS, SUBSTRACTPOLYGONS, INTERSECTPOLYGONS
            PQ = polyclip(p,q,2,varargin{:});
        end
        function r = offsetpolygon(p,offset,varargin)
        % r = OFFSETPOLYGON(p,offset,jointype,par,scale)
		% Alias for POLYOUT, expand/contract a polygon by by a distance 'offset'
        % See also: POLYOUT
        
            r = polyout(p,offset,varargin{:});
        end 
        
        % Other polygon operations
        function Pr = polyrotate(P,R)
		% PolyRotate(P,R) Rotate polygon(s) P using rotation matrix R
		% P can be an array of polygons, with origin at rotation point.
		% R is a 2x2 rotation matrix, or scalar angle (degrees)
		% Result is a rotated polygon object
        
            if isscalar(R), R = [cosd(R),-sind(R);sind(R),cosd(R)]; 
            else, assert(isequal(size(R),[2,2]),'polyrotate:R','Expecting scalar or 2x2 matrix'); end
            
			if numel(P) > 1
				for j = numel(P):-1:1
					Pr(j) = polyrotate(P(j),R);
				end
            else
                if size(P.x,1) > 1, P.x = P.x'; end
                if size(P.y,1) > 1, P.y = P.y'; end
				Pr = polygon(R*[P.x;P.y]);
			end
        end
        function P = polytranslate(P,q)
		% Translate polygon P by 2-vector q.
        
            for j = 1:numel(P)
                P(j).x = P(j).x + q(1);
                P(j).y = P(j).y + q(2);
            end
        end
        function P = clipwithplane(P,n,q)
        % C = CLIPWITHPLANE(P,n,q) clip a (set of) 2D-polygon(s) P with a semi-plane of normal n
        % that goes through point q. The result are any parts of P on the positive-normal side of
        % the plane.

			if nargin < 3, q = [0 0 0]'; end
			n = n(:); q = q(:); % make sure vectors are vertical
            n = n/rssq(n);
			
            r = max(hypot([P.x]-q(1),[P.y]-q(2)));
            mask = polygon([n';n(2),-n(1)]*[0,0;0,1;1,1;1,-1;0,-1]'*2*r+q);
            P = intersectpolygons(P,mask);
        end
        function P = refinepoly(P,d0,dist)
        % P = REFINEPOLY(P,d0) - Insert intermediate points along the edges of the polygon(s) P, 
        %   such that the distance between consecutive vertices is allways shorter than d0.
        %
        % P = REFINEPOLY(P,d0,f) - Provide a custom metric for distance f = @(dx,dy) instead of
        %   the default Euclidean distance f = @hypot
        %
        %   NOTE: whatever the distance metric, it is assumed that inserting n points in any
        %   segment will approximately divide the distance between vertices in that segment by
        %   1/(n+1). Implement recursion externally if this might not be the case!
        %
		% See also: CLEANPOLYGON
        
            narginchk(2,3);
            if nargin < 3, dist = @hypot; end

            for k = 1:numel(P)
                A = [P(k).x;P(k).y];                    % list of vertices
                B = circshift(A,-1,2);                  % B is the following vertex for each A
                d = dist(A(1,:)-B(1,:),A(2,:)-B(2,:));  % edge lengths
                n = numel(d);
                
                if any(d > d0)
                    idx = zeros(1,n);       % will serve to place the original vertices
                    addpts = cell(1,n);     % stores additional vertices
                    addidx = cell(1,n);     % serves to place additional vertices

                    for j = find(d > d0)
                    % For segments where new vertices must be inserted...
                    % Add nn points so that distance between points < d0
                        nn = floor(d(j)/d0);                
                        f = (1:nn)/(nn+1);
                        addpts{j} = A(:,j)*f + B(:,j)*(1-f);
                        idx(j) = nn;
                        addidx{j} = (1:nn)+ j + sum(idx(1:j-1)); % know where to place new vertices
                    end
                    addpts = [addpts{:}];
                    addidx = [addidx{:}];

                    idx(2:n) = (2:n)+cumsum(idx(1:end-1)); idx(1) = 1;  % where to place originals
                    A(:,idx) = A;                                       % relocate original vertices
                    A(:,addidx) = addpts;                               % ... and paste new
                    
                    P(k).x = A(1,:); P(k).y = A(2,:);
                end
            end
        end
        function P = cleanpolygon(P,tol)
        % P = CLEANPOLYGON(P,TOL) - Removes any 'unnecesary' vertices in polygon P, within a 
        % tolerance TOL. This means:
        %   - Remove vertices that line within a distance TOL of a line joining two other vertices.
        %   - Remove any resulting polygons that have less than three vertices
		% See also: REFINEPOLY
		
            if isempty(P), return; end
            
            keepers = arrayfun(@(p) numel(p.x) > 2,P);
            if ~any(keepers), P = polygon.empty; return; end
            
            if nargin < 2, tol = getSimOption('minabstol');
            elseif tol == 0
            % For cero tolerance, return polygon(s) with unique vertices
                for j = find(keepers)
                    V = unique([P(j).x;P(j).y]','rows','stable');
                    if size(V,1) < 3, keepers(j) = false;
                    else, P(j).x = V(:,1)'; P(j).y = V(:,2)';
                    end
                end
                P = P(keepers);
                return
            end
            
            for j = find(keepers)
                n = numel(P(j).x);
                if n < 3, keepers(j) = false; continue; end
                
                xy = [P(j).x([1:end,1])',P(j).y([1:end,1])'];
                xy = mdpwl.douglaspeucker(xy,tol);
                if size(xy,1) > 3
                    P(j).x = xy(1:end-1,1)';
                    P(j).y = xy(1:end-1,2)';
                else
                    keepers(j) = false;
                end
            end
            
            if ~any(keepers), P = polygon.empty(); else, P = P(keepers); end
        end
        function E = polyoutline(P,varargin)
            E = polygon.outline([P.x]',[P.y]',varargin{:});
            E = mergepolygons(E,P);
        end
        
        % Plotting
		function varargout = polyplot(P,frontcolor,edgecolor,varargin)
        % H = POLYPLOT(P,frontcolor,edgecolor)
        % H = POLYPLOT(P,frontcolor,edgecolor,'PropertyName',PropertyValue...)
        % H = POLYPLOT(...,AX)
        %
		% Plot a complex-polygon object (a set of non-intersecting polygons) as patches over gcf().
        %
        % P: polygon or array of non-intersecting polygons (for polygons with holes). 
        %   CAUTION: the function isInterior used to determine what is and what's not a hole
        %   doesn't care about p.hole flags, but follows the even-odd fill-rule. Any output of
        %   polygon operations (e.g. mergepolygons, intersectpolygons, etc.) should do the trick, 
        %   so, when in doubt, first run P = mergepolygons(P,polygon()).
        %
        % frontcolor, edgecolor: can be scalar indices for the current ColorOrder, character keys 
        %   for basic colors {'y','m','c','r','g','b','w','k'} or 3/4 vectors for (R,G,B,[alpha]).
        %   Default is the next color in ax.ColorOrder, with a black edge.
        %
        % property-value pairs: will be passed to patch() after x,y, and colors. 
        %   NOTE: Color information passed as part of the property-value pairs will have priority
        %   over frontcolor & edgecolor. e.g. polyplot(p,'r','k','facecolor','b') will plot a blue
        %   (not red) polygon. However polyplot(p,'r','k','facealpha',0.5) will use both arguments,
        %   and is equivalent to polyplot(p,[1 0 0 0.5]).
        %
        % See also: PATCH, POLYGON, TRIANGULATE
      
            if isvoid(P)
            % nothing to plot
				if nargout > 0, varargout{1} = []; end
				return; 
            end
            if ~isempty(varargin) && ishandle(varargin{end}) && isa(varargin{end},'matlab.graphics.axis.Axes')
                ax = varargin{end};
                varargin(end) = [];
            else
                ax = gca(); 
            end
			
			% Intepret color information
			if nargin<3, edgecolor = 'k'; end
            if nargin<2
            % Get next ColorOrder index, based on the number of 'stuff' already plotted in gca
                frontcolor = mod(numel(get(gca,'Children')),size(get(gca,'ColorOrder'),1))+1; 
            end
 
            % It is also possible that color information is provided in varargin...
            % try to get it, leave empty when missing.
            [opt,remargs] = getpairedoptions(varargin,...
                {'facecolor','facealpha','edgecolor','edgealpha'},{[],[],[],[]});
			[cf,alphaf] = sortcolor(frontcolor,opt.facecolor,opt.facealpha); 
            [ce,alphae] = sortcolor(edgecolor,opt.edgecolor,opt.facealpha);

			colorargs = {'FaceColor',cf,'FaceAlpha',alphaf,'EdgeColor',ce,'EdgeAlpha',alphae};
			
            [V,~,F] = poly2vef(P,true);
            if numel(P) > 1 && any([P.hole])
            % if there are holes, use Constrained Delaunay triangulation to find a set of
            % triangles that represent the actual 'filled' area...

                warning_reseter = naptime({'MATLAB:delaunayTriangulation:ConsConsSplitWarnId',... 
                                           'MATLAB:delaunayTriangulation:DupConsWarnId'}); %#ok<NASGU>
                if ~isa(P,'polygon3d')
                % If there are self-intersections, the number of vertices can change, hence the use
                % a different variable VT.
                     %[VT,~,FT,Issues] = triangulate(P);
                     [VT,~,FT] = triangulate(P);
                else
                % In 3D polygons, perform triangulation over the polygon's plane...
                    [Q,R,s] = flatten(P);
                    [VT,~,FT,Issues] = triangulate(Q);
                    
                    % if the number of vertices did not change, avoid precision issues by keeping V
                    if isempty(Issues.V), VT = V;
                    % otherwise transform the resulting vertices... for truly flat polygons the
                    % result should be very close
                    else, VT(:,3) = 0; VT = bsxfun(@plus,R*VT',s)';
                    end
                end
               
                % Plot the triangles, passing face-color info, with no edges
                H(2) = patch(ax,'Faces',FT,'Vertices',VT,colorargs{1:4},'EdgeColor','none',remargs{:});
                % Then plot the edges without faces...
                H(1) = patch(ax,'Faces',F,'Vertices',V,colorargs{5:8},'FaceColor','none',remargs{:});
                
                % % Plot self-intersection issues
                % washold = ishold(); hold on
                % if ~isempty(Issues.V), plot(Issues.V(:,1),Issues.V(:,2),'ro'); end
                % if ~isempty(Issues.E)
                %     edg = reshape(VT(Issues.E',:)',2,2,[]);
                %     edg(:,3,:) = NaN;
                %     edg = reshape(edg,2,[])';
                %     plot(edg(:,1),edg(:,2),'r-','linewidth',2)
                % end
                % if ~washold, hold off; end
            else
            % With no holes, just plot faces and edges as they are...
                H = patch(ax,'Faces',F,'Vertices',V,colorargs{:},remargs{:});
            end
            
            if nargout > 0, varargout{1} = H; end
            
            function [cdata,alpha] = sortcolor(colorinfo,forcecolor,forcealpha)
            % Guess what the user means by colorinfo
            
                [cdata,alpha] = rgba(colorinfo);
                
                % overwrite non-empty forcecolor/forcealpha
                if nargin > 1 && ~isempty(forcecolor), cdata = forcecolor; end
                if nargin > 2 && ~isempty(forcealpha), alpha = forcealpha; end

                assert(~(isscalar(cdata)&&isnan(cdata)),'polyplot:color','Unrecognized color format');
                assert(size(cdata,1)==1,'polyplot:color','polygon:polyplot requires scalar colors');
            end
        end
        function [h,U] = labelpolygons(P,lbl,varargin)
        % [H,U] = LABELPOLYGONS(P,LBL,varargin) - place labels LBL (N cell-str) on N array of
        %   (simple) polygons P, trying to avoid overlap among labels.
        %
        %   Additional arguments are passed to TEXT().
        %   Returns TEXT handles H, and optionally an array U of rectangles with label extents.
        %
        % TODO: replace by an external, standardized placement library?
        % TODO: handle complex polygons and 3D polygons
        % TODO: use distance to straight-skeleton instead of random polyon points?
        % TODO: turn hard-coded parameters into optional inputs

            NPTS = 32;     % number of random points in polygon that will 'attract' labels
            EXP_X = 4;     % decay exponent for 'repulsive force' among labels
            EXP_C = 4;     % decay exponent for 'attractive force' from labels to polygon
            MIN_R2X = 1/4; % min. square distance, label to label (clips 1/r^EXP_X)
            MIN_R2C = 1/4; % min. square distance, label to polygon (clips 1/r^EXP_C)
            MAX_MOV = 1/8; % max. step displacement, in units of label extents
            
            MAX_ITER = 10;
            TOL = 1e-3;

            N = numel(P);
            assert(numel(lbl) == N);

            [Q,C,rc] = arrayfun(@(p) random_points(p,NPTS),P,'UniformOutput',false);

            Q = permute(cat(3,Q{:}),[3,2,1]);   % random 'attraction points' inside each polygon
            C = cat(2,C{:})';                   % centroids
            rc = cat(1,rc{:});                  % characteristic lenghts

            % plot labels at centroids, just to measure their extents
            h = text(C(:,1),C(:,2),zeros(N,1),lbl,varargin{:});
            ex = cat(1,h.Extent);
            rx = ex(:,3)/2;
            ry = ex(:,4)/2;
            offset = ex(:,1:2) + [rx,ry] - C;
            X = C;

            % optimize label positions
            S = permute(cat(3,1./(rx + rx'),1./(ry + ry')),[1,3,2]);
            for l = 1:MAX_ITER
                Dx = (X - permute(X,[3,2,1])).*S;
                r2x = max(MIN_R2X,dot(Dx,Dx,2));
                Dc = (X - Q)./rc;
                r2c = max(MIN_R2C,dot(Dc,Dc,2));

                f = sum(Dx.*r2x.^(-EXP_X/2),3) - mean(Dc.*r2c.^(-EXP_C/2),3);
                f = f.*rc;
                e = sign(f).*min(abs(f),[rx,ry].*MAX_MOV);

                if all(abs(e) < [rx,ry].*TOL,'all'), break; end
                X = X + e;
            end

            % move labels to their optimized position
            for j = 1:N, set(h(j),'Position',[X(j,:)-offset(j,:),0.1]); end

            if nargout > 1
                U = arrayfun(@(j) polytranslate(polygon(2*rx(j),2*ry(j)),X(j,:)),1:N);
            end
        end
        function [pts,c,r] = random_points(p,N,THRESH)
        % [PTS,C,R] = RANDOM_POINTS(P,N) - return N random points inside polygon p or,if it is  
        % closely degenerate, return random points from its outline.
        %
        % TODO: distribute points using delaunay triangulation and barycentric coordinates?
        
            if nargin < 3 || isempty(THRESH), THRESH = 1e-2; end

            L = p.perimeter;
            c = centroid(p);
            r = 2*p.area/L;

            if r/L < THRESH
                pts = zeros(0,2);
            else
                [w,h,c] = rectangleproperties(p);
                q = cell(100,1);
                n = 0;
                for j = 1:100
                    pts = (2*rand(N,2)-1).*[w,h]+c';
                    in = insidepolygon(p,pts(:,1),pts(:,2));
                    q{j} = pts(in,:);
                    n = n + size(q{j},1);
                    if n > N, break; end
                end
                pts = cat(1,q{:});
                pts(N+1:end,:) = [];
            end

            if size(pts,1) < N
                p = refinepoly(p,r/2);
                V = poly2vef(p);
                pts = [pts;V(randperm(size(V,1)),:)];
                pts(N+1:end,:) = [];
            end

            for j = size(pts,1)+1:N, pts(j,:) = c; end
        end
    end
    
    methods (Static)
    % Unpack is static as it currently takes a structure
    
        function p = unpack(S)
        % S = UNPACK(P) convert polygon P CLIPPER ouput coordinates from int62 values
        % UNPACK is called by default inside POLYCLIP, unless input polygons are pre-packed
        % See also: POLYCLIP, CLIPPER, PACK
        
            if isempty(S), p = polygon.empty; return; end
            
            % if nargin < 2 || isempty(scale)
            %     if isfield(S,'scale')
            %         scale = double([S.scale]);
            %     else
            %         switch(class(S.x))
            %         case 'int16'
            %             scale = polygon.SCALE16;
            %         case 'int64'
            %             scale = polygon.SCALE16;
            %         otherwise
            %             error('No default scale for class: %s',class(P.x));
            %         end
            %     end
            % else
            %     scale = double(scale);
            %     if isfield(S,'scale')
            %         scale = double([S.scale]);
            %     end
            % end
            % if isscalar(scale) && numel(S) > 1, scale = repmat(scale,size(S)); end
            scale = double([S.scale]);

            p(numel(S)) = polygon();
            for j = 1:numel(S)
                p(j).x = double(S(j).x)/scale(j);
                p(j).y = double(S(j).y)/scale(j);
                % n = numel(p(j).x);
                % Equivalent to fixorientation(polygon(S(j)),true):
                % p(j).hole = sum(p(j).x.*p(j).y([2:n,1]) - p(j).x([2:n,1]).*p(j).y) < 0;
                p(j).hole = [];
                p(j).hole = p(j).signedarea < 0;
            end
        end
        function in = inpackedpolygon(P,x,y)
        % IN = INPACKEDPOLYGON(P,X,Y) - equivalent to INSIDEPOLYGON(UNPACK(P),X,Y)
        
            if isempty(P), in = false(size(x)); return; end
            
            % if nargin < 4 || isempty(scale)
            %     switch(class(P.x))
            %     case 'int16'
            %         scale = polygon.SCALE16;
            %     case 'int64'
            %         scale = polygon.SCALE16;
            %     otherwise
            %         error('No default scale for class: %s',class(P.x));
            %     end
            % end
            scale = double([P.scale]);

			for j = 1:numel(P)
                px = double(P(j).x)/scale(j);
                py = double(P(j).y)/scale(j);
                
                if ~isfield(P(j),'hole')
                    xi = px - mean(px);
                    yi = py - mean(py);
                    xj = circshift(xi,-1);
                    yj = circshift(yi,-1); 
                    P(j).hole = sum((xi-xj).*yj - (yi-yj).*xj) < 0;
                end
				if P(j).hole
					if j == 1, in = true(size(x)); end
					in = in & ~inpolygon(x,y,px,py);
				else
					if j == 1, in = false(size(x)); end
					in = in | inpolygon(x,y,px,py);
				end
			end
        end
        function A = packedarea(P)
		% A = PACKEDAREA(P) - equivalent to POLYGONAREA(UNPACK(P))
        
            A = 0;
            if isempty(P), return; end
            
            % if nargin < 2 || isempty(scale)
            %     switch(class(P.x))
            %     case 'int16'
            %         scale = polygon.SCALE16;
            %     case 'int64'
            %         scale = polygon.SCALE16;
            %     otherwise
            %         error('No default scale for class: %s',class(P.x));
            %     end
            % end
			
            %k = 2*double(scale)^2;
            k = 2*double([P.scale]).^2;
            for j = 1:numel(P)
                if isempty(P(j)) || numel(P(j).x) < 3, continue; end
                xi = double(P(j).x); xi = xi - mean(xi);
                yi = double(P(j).y); yi = yi - mean(yi);
                xj = circshift(xi,-1);  % faster than x([2:end,1])
                yj = circshift(yi,-1); 
                % A = sum(xi·yj - yj·xi)/2 = sum[(xi-xj)·yj - (yi-yj)·xj]/2, better precision
                A = A + sum((xi-xj).*yj - (yi-yj).*xj)/k(j);
                % A = A + sum(xi.*yi([2:end,1]) - xi([2:end,1]).*yi)/k;
            end
        end
        function [P,T,a] = outline(x0,y0,offset,n)
        % P = OUTLINE(X,Y) - Return a non-convex polygonal boundary P of points (X,Y). Namely,
        %   P is the outline of the ALPHASHAPE of X,Y with a radius ALPHA such that P is a simple
        %   polygon (no holes, unique vertices, non-self-intersecting, and containing all points
        %   X,Y). Similar to BOUNDARY, but with an adaptive radius that ensures a simple polygon.
        %   
        % P = OUTLINE(X,Y,OFFSET,N) - To guarantee that the resulting polyong has non-zero
        %   thickness everywhere (e.g. colinear points), an offset can be applied to all points, 
        %   forming an N-sided polygon of radius OFFSET around each original point, before finding
        %   the alpha-shape. OFFSET is also set as the minimum alpha-radius.
        %
        % [P,T,ALPHA] = OUTLINE(...) - return additionally the ALPHATRIANGULATION and alpha-radius
        %
        % NOTE: for a POLYGON Q, P = OUTLINE([Q.x],[Q.y],...) does not ensure that P will fully
        %   contain Q. Use POLYOUTLINE(Q) for that purpose, which merges the result of OUTLINE with
        %   the original POLYGON.
        %
        % See also: POLYGON, ALPHASHAPE, CONVEXHULL, POLYOUTLINE, BOUNDARY

            if nargin < 3 || isempty(offset), offset = 0; end
            if nargin < 4 || isempty(n), n = 6; end
            isnice = @(x) isscalar(x) && isfinite(x) && isreal(x) && x >= 0;
            assert(isnice(offset),'Bad offset');
            assert(isnice(n) && mod(n,1) == 0,'Bad n');

            if offset > 0
               p = polygon(n);
               x = x0(:) + offset*[0,p.x];
               y = y0(:) + offset*[0,p.y];
            else
               x = x0(:);
               y = y0(:);
            end
            [P,idx] = unique([x(:),y(:)],'rows','stable');
            if ~isequal(idx',1:numel(x)), x = x(idx); y = y(idx); end
            
            assert(numel(x) > 2 && rank(P) >= 2,'rank([x,y] < 2, non-zero offset required');

            shp = alphaShape(double(P),Inf);
            shp.HoleThreshold = shp.area();
            if offset > 0
            % Set alpha to the next-bigger-than-OFFSET value in alphaSpectrum 
                aa = alphaSpectrum(shp);
                a = aa(find(aa > offset,1,'last'));
                if isempty(a), a = offset; end  % offset > max(aa), should yield convex-hull
            else
            % Start with the lowest alpha value for which the boundary is a single region
                a = shp.criticalAlpha('one-region');
            end

            if ~yieldsniceboundary(a)
            % search within critical alpha and max. alpha (convex-hull)
                aa = alphaSpectrum(shp);
                fun = @(j) yieldsniceboundary(aa(round(j)))*2-1;
                j = bisection(fun,1,find(aa == a),[0.8,0]);
                j = round(j-0.5);
                a = aa(j);
            end
            
            if ~yieldsniceboundary(a)
                warning('Failed to find non-convex outline, returning convex-hull');
                b = convhull(x,y,'simplify',true);
            else
                b = shp.boundaryFacets();
            end
            P = polygon(x(b(:,1)),y(b(:,1)));
            
            if nargout > 1
                T = triangulation(shp.alphaTriangulation(),shp.Points);
            end
            
            function isok = yieldsniceboundary(a)
            % set alhpa = a, and check if boundary is a simple polygon
            
                shp.Alpha = a;
                b = shp.boundaryFacets();

                isok = all(b(:,1) == circshift(b(:,2),1)) && ...     % single loop
                       all(nonzeros(accumarray(b(:),1)) == 2) && ... % non-repeated vertices
                       all(inpolygon(x0(:),y0(:),x(b(:,1)),y(b(:,1))));
            end
        end

    end
end
