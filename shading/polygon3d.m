 classdef polygon3d < polygon
 % POLYGON3D class - used to hold both rotated planar polygons, and generic sets of
 %  3D vertices (e.g. spherical polygons).
 %
 % TODO: POLYGON should be a subclass of POLYGON3D, and not the other way around

    properties
		z
	end
    methods
        function obj = polygon3d(varargin)
		% POLYGON3D()	- Empty polygon
		% POLYGON3D(S) - Polygon from [array of] structure(s) or 2D polygon(s).
		%				 S(j).x, S(j).y are vectors with the vertices of a poligon j.
		%				 Holes (cut-out regions) are ignored.
		% POLYGON3D(XYZ) - Polygon from 3xn array, where every XYZ(:,j) is a vertex.
		% POLYGON3D(X,Y,Z) - 3D polygon, with vertices at X,Y,Z.
		% POLYGON3D(az,el,r,'sph') - 3D polgon from spherical coordinates (angles in degrees)
		% POLYGON3D(th,r,z,'pol') - 3D polygon from cylindrical coordinates (angles in degrees)
        % NOTE: for all coordinate systems, all arrays must have the same number of elements or
        %   or be scalar (implying a constant dimension). Examples:
        %   POLYGON3D([1 1 -1 -1],0,[1 -1 -1 1]) is a square on the XZ plane
        %   POLYGON3D(0:60:300,1,2,'pol') is a unit hexagon on the plane Z = 2.
		
            narginchk(0,4);
            argerr = {'polygon3d:construct','Cannot recognize arguments for constructor'};
            switch(nargin)
            case 0
            % empty polygon constructor 
                obj.x = [];
                obj.y = [];
                obj.z = [];
                obj.hole = false(0,0);
            case 1
                P = varargin{1};
                if isnumeric(P) 
                % from set of points [Px;Py;Pz]
                    if size(P,1) ~= 3, P = P'; end
                    assert(size(P,1) == 3,argerr{:});
                    obj.x  = P(1,:); obj.y = P(2,:); obj.z = P(3,:); obj.hole = false;
                else
                % 3D polygon from structure or 2D polygon (or array of any)
                    assert(isa(P,'polygon') || isstruct(P) && all(isfield(P,{'x','y','hole'})),argerr{:});
                    if isempty(P), P = polygon3d.empty; end
                    
                    s = num2cell(size(P));
                    obj(s{:}) = polygon3d();
                    for j = 1:numel(P)
                        obj(j).x = P(j).x;
                        obj(j).y = P(j).y;
                        obj(j).z = zeros(size(P(j).x));
                        obj(j).hole = P(j).hole;
                    end
                end
            case {3,4}
                
                % Get coordinate system
                if nargin < 4, varargin{4} = 'car';
                else, assert(ischar(varargin{4}),argerr{:});
                end

                % Convert to cartesian coordinates, if required
                switch lower(varargin{4})
                case 'sph'
                    [varargin{1:3}] = sph2cart(varargin{1}*pi/180,varargin{2}*pi/180,varargin{3});
                case 'pol'
                    [varargin{1:3}] = pol2cart(varargin{1}*pi/180,varargin{2},varargin{3});
                case 'car'
                    % do nothing
                otherwise
                    error(argerr{:});
                end

                % Check that all dimensions have the same size. Expand scalars, if required
                n = max(cellfun(@numel,varargin(1:3)));
                for j = 1:3
                    switch numel(varargin{j})
                    case n, varargin{j} = varargin{j}(:)';
                    case 1, varargin{j} = ones(1,n)*varargin{j};
                    otherwise
                        error('polygon3d:arrsize','All arrays must be equal (or scalar)');
                    end
                end

                obj.x = varargin{1};
                obj.y = varargin{2};
                obj.z = varargin{3}; 
                obj.hole = false;
            end
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
                    S(j).z = P(j).z;
                    S(j).hole = P(j).hole;
                end
            end
        end
      	function [V,E,F] = poly2vef(p,unif)
        % [V,E,F] = POLY2VEF(P,UNIF) - Generates a list of vertices, edges, and faces of a (set of) 
        % polygon(s) P.
        %
        %   V: Nv·3 array of vertices, where V(j,:) = [xj,yj]. When numel(p) > 1, the vertices
        %      of all polygons p(k) will be concatenated in order.
        %   E: Ne·2 array of edges (vertex-index-pairs), where E(k,:) = [i,j] implies that there is 
        %      a vertex from V(i,:) to V(j,:).
        %   F: When UNIF = false (default), F is a cell-array of vertex-index-tuples. Each element 
        %      k of F represents a face formed by vertices F{k} = [a,b,..]
        %      When UNIF = true, F is packed into a uniform Nf·t array of (vertex-index-tuples), 
        %      where t is the maximum number of vertices for all faces, and face k is given by 
        %      F(k,~isnan(F(k,:)) i.e. for faces with less than t vertices, non-used positions are
        %      filled with NaN.
        %
        % See also: POLYGON3D.VF2POLY
        
            if nargin < 2, unif = false; end
        
            Np = numel(p);
            V = cell(Np,1);
            Nv = zeros(Np,1);
            for j = 1:Np
                V{j} = [p(j).x;p(j).y;p(j).z]';
                Nv(j) = numel(p(j).x);
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
        
        function varargout = flatpolymethod(P,fun,nout,varargin)
            if nargin < 3, nout = 1; end
            [Q,R,c] = flatten(P);
            for j = 1:numel(varargin)
                if isa(varargin{j},'polygon3d')
                    varargin{j} = flatten(varargin{j},R,c);
                end
            end
            [varargout{1:nout}] = fun(Q,varargin{:});
            for j = 1:nout
                if isa(varargout{j},'polygon')
                    Q = polygon3d(varargout{j});
                    varargout{j} = polyrotate(Q,R);
                end
            end
        end
        
        function c = centroid(p)
            q = polygon();
            if all(p.x == p.x(1))
                q.x = p.y;
                q.y = p.z;
                c = [p.x(1);centroid(q)];
            elseif all(p.y == p.y(1))
                q.x = p.x;
                q.y = p.z;
                cxz = centroid(q);
                c = [cxz(1);p.y(1);cxz(2)];
            elseif all(p.z == p.z(1))
                q.x = p.x;
                q.y = p.y;
                c = [centroid(q);p.z(1)];
            else
                q.x = p.x;
                q.y = p.z;
                cxz = centroid(q);
                q.y = p.y;
                c = [centroid(q);cxz(2)];
            end
        end
        
        % Other polygon operations
        function P = polyrotate(P,varargin)
		% Q = POLYROTATE(P,R) - Rotate polygon(s) P using 3x3 rotation matrix R (left-multiplied)
        % Q = POLYROTATE(P,ANG,CONV) - Rotate by (degree) angles ANG = [A,B,C..] around the 
        %    (rotated) axes specified in CONV (e.g. 'ZXZ'). See ROTMAT for details.
        % Q = POLYROTATE(..,S) - Rotate around point S (default is origin [0;0;0]);
        %
		% See also: ROTMAT
        
            narginchk(2,4);
            if isempty(P), return; end
            
            S = varargin{end};
            if isnumeric(S) && numel(S)==3  % S is point of rotation
                S = S(:); 
                varargin = varargin(1:end-1);
            else, S = []; 
            end
            
            switch numel(varargin)
                case 0, R = []; % error
                case 1
                    R = varargin{1};
                    if ~isnumeric(R) || ~isequal(size(R),[3,3]), R = []; end  % error
                otherwise
                    try R = polygon3d.rotmat(varargin{:}); catch, R = []; end
                    assert(size(R,3) == 1,'POLYROTATE does not (yet) work for multiple-rotations');
            end
            assert(~isempty(R),'polyrotate:args','Arguments do not match any recognized pattern.');
            
            for j = 1:numel(P)
                V = [P(j).x;P(j).y;P(j).z];
                if ~isempty(S), V = bsxfun(@minus,V,S); end
                V = R*V;
                if ~isempty(S), V = bsxfun(@plus,V,S); end
                P(j).x = V(1,:); P(j).y = V(2,:); P(j).z = V(3,:);
            end
        end
		function P = polytranslate(P,q)
		% Translate polygon P by 3-vector q.
            for j = 1:numel(P)
                P(j).x = P(j).x + q(1);
                P(j).y = P(j).y + q(2);
                P(j).z = P(j).z + q(3);
            end
		end
        function P = clipwithplane(P,n,q)
        % C = CLIPWITHPLANE(P,n,q) clips a (set of) 3D-polygon(s) P using a plane with normal n
        % that goes through point q. The result are any parts of P on the positive-normal side of
        % the plane.
        
			if nargin < 3, q = [0 0 0]'; end
			n = n(:); q = q(:); % make sure vectors are vertical
			
            preserve = true(size(P));
            for k = 1:numel(P)
                V = [P(k).x;P(k).y;P(k).z]; 
                zz = n'*bsxfun(@minus,V,q); % magnitude of normal distances from cutting plane
                ap = (zz >= 0);             % above-plane flags
                
                % There are three possibilities:
                % 1. fully in front of plane, do nothing to P(k)
                    if all(ap), continue; end
                    
                % 2. fully behind plane (polygon will be deleted)
                    if ~any(ap), preserve(k) = false; continue; end 
                
                % 3. polygon is actually clipped...
                npt = size(V,2);
                cuts = xor(ap,ap([2:end,1])); % get index of edges that cross the plane
                
                % Make an index of positions for the original vertices, leaving gaps for new ones
                % e.g. if cuts = [0,0,1,0...], idx = [1,2,3,5,6...]
                idx = ones(1,npt);
                idx(2:npt) = (2:npt)+cumsum(cuts(1:end-1));
                
                % Calculate intersections, append to the list of vertices...
                a = find(cuts); b = mod(a,npt)+1; nc = numel(a);
                V(:,end+1:end+nc) = V(:,a)+bsxfun(@times,-zz(a)./(zz(b)-zz(a)),V(:,b)-V(:,a));
                
                % ... and update index (point new vertices to the gaps left originally)
                idx(end+1:end+nc) = idx(a)+1;
                idx(1:npt) = idx(1:npt).*ap;  % leave zero for points below plane
                
                % Place nonzero-vertices where they should be, and remove the rest
                V(:,nonzeros(idx)) = V(:,idx > 0);
                V = V(:,sort(nonzeros(idx)));

                P(k).x = V(1,:); P(k).y = V(2,:); P(k).z = V(3,:);
            end
            P = P(preserve);
        end
        function P = refinepoly(P,d0)
        % P = REFINEPOLY(P,d0) - Insert intermediate points along the edges of the polygon(s) P, 
        % such that the distance between consecutive vertices is allways shorter than d0.
		% See also: CLEANPOLYGON

            for k = 1:numel(P)
                A = [P(k).x;P(k).y;P(k).z];        % get list of vertices
                B = A(:,[2:end,1]);                % B is the following vertex for each A
                d = rssq(A-B,1);                   % edge lengths
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
                    
                    P(k).x = A(1,:); P(k).y = A(2,:); P(k).z = A(3,:);
                end
            end
        end
        function P = cleanpolygon(P,tol)
		% P = CLEANPOLYGON(P,TOL) - Removes any 'unnecesary' vertices in polygon P, within a 
        % tolerance TOL. This means:
        %   - Remove vertices that lie within a distance TOL of a line joining two other vertices.
        %   - Remove any resulting polygons that have less than three vertices
		% See also: REFINEPOLY
            
            if isempty(P), return; end
            
            keepers = arrayfun(@(p) numel(p.x) > 2,P);
            if ~any(keepers), P = polygon3d.empty; return; end
            
            if nargin < 2, tol = getSimOption('minabstol');
            elseif tol == 0
            % For cero tolerance, return polygon(s) with unique vertices
                for j = find(keepers)
                    V = unique([P(j).x;P(j).y;P(j).z]','rows','stable');
                    if size(V,1) < 3, keepers(j) = false;
                    else, P(j).x = V(:,1)'; P(j).y = V(:,2)'; P(j).z = V(:,3)';
                    end
                end
                P = P(keepers);
                return
            end

            for j = find(keepers)
                n = numel(P.x);
                if n < 3, keepers(j) = false; continue; end

                V = [P(j).x;P(j).y;P(j).z];
                useful = true(1,n);
                one2n  = @(x) mod(x-1,n)+1; % used to interprete point n+1 as 1 (closed polygon)

                a = 1; % last pivot point
                b = a+2; % next pivot point 
                while a < n %(*)
                    % make a list of indices for points between a and b (originally just one)
                    jj = one2n(a+1:b-1); 

                    % if all points jj fall close enough to the line a -> b
                    if all(point2line(V(:,a),V(:,one2n(b)),V(:,jj))<tol)
                        % try to go even further, to b+1
                        b = b+1; while ~useful(one2n(b)), b=b+1; end
                        useful(jj) = false;
                        if nnz(useful) < 3, keepers(j) = false; continue; end
                    else
                        % otherwise go back to b-1, and start fresh with a new interval
                        a = b-1;
                        b = a+2; while ~useful(one2n(b)), b=b+1; end
                    end
                end
                P(j).x = P(j).x(useful); P(j).y = P(j).y(useful); P(j).z = P(j).z(useful);
            end
            
            if ~any(keepers), P = polygon3d.empty; else, P = P(keepers); end
            
            function [d,Q] = point2line(a,b,C)
            % Calculates the minimum distance from the points in C to the line a->b
            % a and b must be 3-vectors, C can be a 3·N array

                A = repmat(a,1,size(C,2));
                if all(a==b), d = rssq(C-A); Q = A; return; end % point-to-point distance
                
                u = (b-a)/norm(b-a); % unit vector in the direction a -> b
                Q = (C-A)-u*(u'*(C-A));
                d = rssq(Q);
                Q = Q + C;
            end
		end
      
        function P = polyprojectshadow(P,varargin)
		% Q = POLYPROJECTSHADOW(P,s,n,q) - POLYGON3D.PROJECTSHADOW applied to polygon(s) P
        % Q = POLYPROJECTSHADOW(P,s,T,V)
			
            narginchk(2,4);
            if isempty(P), return; end
            [Vp,~,F] = poly2vef(P);
            Vp = polygon3d.projectshadow(Vp,varargin{:});
            P = polygon3d.vf2poly(Vp,F);
        end
        
        function [n,R,s] = normal(P)
        % Calculate polygon normal, using Principal-Component-Analysis of all vertices in p(:)
        % optionally return the rotation matrix R, and average of vertices s
        
            V = poly2vef(P)';
            s = mean(V,2);
        
            % return NaN for degenerate polygons
            if size(V,2) < 3 || rank(V) < 2, n = nan(3,1); R = nan(3); return; end

            % The idea is to get the least-eigenvector of the covariance matrix XX', i.e. the
            % direction in which the data shows minimal variance.
            % According to (*), using Singular-Value-Decomposition is numerically better than
            % calculating XX' directly: XX'=(USV')(USV')' = US²U'
            % (*) https://math.stackexchange.com/questions/3869/what-is-the-intuitive-
            %     relationship-between-svd-and-pca/3871#3871
            
            V = V-s;           % Substract the average of vertices
            [R,~,~] = svd(V);  % left-singular-vectors (decreasing singular-values)
            n = R(:,3);        % the plane normal is the left-singular-vector for min(diag(S))
        end
        
        function [Q,R,c] = flatten(P,R,c,tol)
        % [Q,R,c] = FLATTEN(P,[R,S,TOL]) - return planar polygon Q = (P-c)R, assuming a rotation 
        %   matrix R and rotation center s can be found (or are provided) such that:
        %   |z| = |(P-c)R·k| < TOL, where k is the unit z vector [0;0;1].
        
            if nargin < 4, tol = eps(1); end
            if nargin < 3, [~,R,c] = normal(P); end
            
            [V,~,F] = poly2vef(P);
            
            V = (V-c')*R;
            if max(abs(V(:,3))) > max(tol,max(eps(V(:))))
                warning('flatten:zz','polygon is not flat, by up to ±%e',max(abs(V(:,3))));
            end
            
            Q = polygon3d.vf2poly(V(:,1:2),F);
            [Q.hole] = deal(P.hole);
        end
	end
	
	methods(Static)
		function [R,ZXZ] = rotmat(ang,CON)
        % R = ROTMAT([A,B,C,D,...],CONV) - returns the resulting rotation matrix for an
		% unlimited series of successive intrinsic rotations. The (rotated) axes around which
		% rotations by (degree) angles A,B,C,... are performed are provided in CONV, as a string 
        % of {'X','Y','Z'} characters, one for each angle.
        %
        % [A,B,C,...] must be a row-vector of m angles, where m = numel(CON), or an n×m array 
        % representing n independent rotations, for which ROTMAT returns a 3×3×n stack of matrices.
        % NOTE: the resulting matri(x/ces) are to be left-multiplied by the points to be rotated.
        %
        % Example: ROTMAT([A,B,C,D,...],'ZXZY...'), performs succesive rotations
		%   A degrees around Z 
		%   B degrees around X' (i.e. X after rotation by A)
		%   C degrees around Z" (i.e. Z after rotation by A and B)
        %   D degrees around Y"' 
        %   ...
        %
        % [~,ZXZ] = ROTMAT([A,B,C,D,...],CONV) - additionally returns an nx3 list of Euler angles
        %   such that ROTMAT([A,B,C,D,...],CONV) == ROTMAT(ZXZ,'ZXZ') within eps(2*pi)
        %
        % See also: POLYROTATE
            
            narginchk(2,2);
            assert(isnumeric(ang) && ismatrix(ang),'First argument must be an array of angles');
            assert(ischar(CON) && numel(CON) == size(ang,2),'CONV must be a string matching the number of angle-columns');
            CON = upper(CON);
            assert(all(CON >= 'X' & CON <= 'Z'),'CONV string contains something other than X,Y,Z');

            ca = permute(cosd(ang),[2 3 1]);    % [m,1,n] stack columns of cosines
            sa = permute(sind(ang),[2 3 1]);    % [m,1,n] stack coumns of sines
            C = @(j) ca(j,1,:);                 % C(j) = [1,1,n] stack of cos(ang(:,j))
            S = @(j) sa(j,1,:);                 % ...
            
            switch CON
            % Common cases, premultiplied for speed
                case 'ZXZ'    
                    R = [C(1).*C(3)-C(2).*S(1).*S(3), - C(1).*S(3)-C(2).*C(3).*S(1),  S(1).*S(2);
                         C(3).*S(1)+C(1).*C(2).*S(3),  C(1).*C(2).*C(3)-S(1).*S(3) , -C(1).*S(2);
                                 S(2).*S(3)         ,            C(3).*S(2)        ,    C(2)   ];
                case 'ZYX'
                    R = [C(1).*C(2), C(1).*S(2).*S(3)-C(3).*S(1),  S(1).*S(3)+C(1).*C(3).*S(2);
                         C(2).*S(1), C(1).*C(3)+S(1).*S(2).*S(3),  C(3).*S(1).*S(2)-C(1).*S(3);
                            -S(2)  ,           C(2).*S(3)         ,           C(2).*C(3)     ];
                        
                % FUTURE: Wikipedia has the complete list...
                % https://en.wikipedia.org/wiki/Euler_angles#Intrinsic_rotations
                
                otherwise
                % for any other convention, compose by individual rotations
                    R = simplerotmat(C(1),S(1),CON(1));
                    for j = 2:numel(CON)
                        Rj = simplerotmat(C(j),S(j),CON(j));
                        for k = 1:size(ang,1), R(:,:,k) = R(:,:,k)*Rj(:,:,k); end
                    end
            end
            
            if nargout > 1
            % Return true Euler angles for ZXZ intrinsic convention
               ZXZ(:,1) = atan2d(R(1,3,:),-R(2,3,:));
               ZXZ(:,2) = atan2d(hypot(R(1,3,:),R(2,3,:)),R(3,3,:));
               ZXZ(:,3) = atan2d(R(3,1,:),R(3,2,:));
            end
			
			function R = simplerotmat(cosA,sinA,axis)
                R = repmat(eye(3),1,1,numel(cosA));
				switch axis
				case 'X' % R = [1,0,0;0,C,-S;0,S,C];
                    R(2,2,:) = cosA; R(2,3,:) = -sinA; R(3,2,:) = sinA;  R(3,3,:) = cosA;
				case 'Y' % R = [C,0,S;0,1,0;-S,0,C];
                    R(1,1,:) = cosA; R(1,3,:) = sinA;  R(3,1,:) = -sinA; R(3,3,:) = cosA;
				case 'Z' % R = [C,-S,0;S,C,0;0,0,1];
                    R(1,1,:) = cosA; R(1,2,:) = -sinA; R(2,1,:) = sinA;  R(2,2,:) = cosA;
				end
			end
        end
        
        function P = vf2poly(V,F)
        % P = VF2POLY(V,F) - generate polygon-objects from a set of vertices (V) and faces (F).
        %   V: Nv·3 array of vertices, where V(j,:) = [xj,yj,zj].
        %   F: Nf·t array or Nf-cell-array of face indices
        %   P: polygon or array of polygons 
        %      NOTE: all hole-flags will be set to zero!, use P = fixorientation(P,true) to 
        %      asign flags depending on vertex-order.
        %
        % TODO: switch to quaternion representation
        %
        % See also: POLY2VEF
            
            narginchk(2,2);
            is3d = size(V,2) == 3;
            assert(isnumeric(V) && (is3d ||(size(V,2) == 2)),'vf2poly:V',...
                'expecting n x 2 or n x 3 array of vertices as first argument');
            
            if isempty(F), if is3d, P = polygon3d.empty; else, P = polygon.empty; end; end
            
            if is3d, P(size(F,1),1) = polygon3d();
            else, P(size(F,1),1) = polygon();
            end
            
            if iscell(F)
                f = [F{:}];
                if any(mod(f,1)~=0) || any(f > size(V,1)) || any(f < 1) 
                    error('vf2poly:fidx','Index values in F must be in [1:rows(V)] or be NaN');
                end

                for j = 1:numel(F)
                    P(j).x = V(F{j},1)';
                    P(j).y = V(F{j},2)';
                end
                if is3d, for j = 1:numel(F), P(j).z = V(F{j},3)'; end; end
            else
                f = ~isnan(F);
                if any(mod(F(f),1)~=0) || any(F(f) > size(V,1)) || any(F(f) < 1) 
                    error('vf2poly:fidx','Index values in F must be in [1:rows(V)] or be NaN');
                end

                for j = 1:size(F,1)
                    P(j).x = V(F(j,f(j,:)),1)';
                    P(j).y = V(F(j,f(j,:)),2)';
                end
                if is3d, for j = 1:size(F,1), P(j).z = V(F(j,f(j,:)),3)'; end; end
            end
            [P.hole] = deal(0);
        end
        
        function Vp = projectshadow(Vp,s,varargin)
		% Q = PROJECTSHADOW(V,s,n,q), where n and q are 3-vectors, projects vertices V in 
        %   direction s over a plane with normal n that goes trhough point q. The result is 
        %   the shadow/antishadow[!] cast by V over the plane, i.e.
        %
        %       P' = V' - [n·(V'-q)/(n·s)]·s
		%
        %   [!] NOTE: PROJECTSHADOW doesn't check if points actually lie above the plane
        %   i.e. (P-q)·n > 0, or if the points are actually projected towards s and not away 
        %   from s, i.e. [P(j)-Q(j)]·s < 0
        %
        %	s, n, q must be 3-vectors
		%	s - defines the direction of the projection (sign and |s|~=1 don't matter!)
		%   n - defines the normal to the plane (sign and |n|~=1 don't matter!)
		%   q - is a vector from the origin to any point over the plane. If not provided it will
		%		be assumed that a = 0, i.e. the plane passes through the origin.
		%
		% Q = PROJECTSHADOW(P,s,Tri,V), where Tri is a triangulation and V a set of vertices:
		%	Projects polygon P in direction s over the triangular mesh defined by Tri and V.
			
			% input parsing (resolve function overloading)
            narginchk(2,4);
            if isempty(Vp), return; end
            
            validateattributes(s,{'numeric'},{'real','finite','vector','numel',3},'','s');
            s = s(:);

            switch nargin
            case 2, n = [0;0;1]; q = [0;0;0];
            case 3
                if isnumeric(varargin{1})
                    n = varargin{1};
                    q = [0;0;0];
                    projecttoplane = true;
                else
                    validateattributes(varargin{1},{'triangulation'},{'nonempty'},'','T');
                    V = varargin{1}.Points;
                    T = varargin{1}.ConnectivityList;
                    projecttoplane = false;
                end
            case 4
                if numel(varargin{1})==3
                    n = varargin{1}(:);
                    q = varargin{2}(:);
                    projecttoplane = true;
                else
                    T = varargin{1};
                    V = varargin{2};
                    projecttoplane = false;
                end
            end
            
			% project to plane
			if projecttoplane     
                
                validateattributes(n,{'numeric'},{'real'},'','n');
                validateattributes(q,{'numeric'},{'real','numel',3},'','n');
                        
                % Get projected verices: Q = P - [n·(P-q)/(n·s)]·s (*)
                Vp = Vp-((Vp-q')*n/dot(n,s))*s';
                return;
			end
			
			% project to mesh
            validateattributes(V,{'numeric'},{'real','2d','size',[NaN,3]},'','V');
            validateattributes(T,{'numeric'},{'2d','integer','positive','size',[NaN,3]},'','T');

			options.triangles = 'one sided';	% only intersections in single direction
			options.ray = 'ray'; 
			options.border = 'inclusive';
			options.epsilon = 1e-5;
			
			nT = size(T,1);
            nP = size(Vp,1);
            
            nchunks = ceil((nP*nT*64)/maxarraysize());
            nJ = ceil(nP/nchunks);
            
            v1 = single(repelem(V(T(:,1),:),nJ,1));
            v2 = single(repelem(V(T(:,2),:),nJ,1));
            v3 = single(repelem(V(T(:,3),:),nJ,1));
            ss = single(repmat(s',nT*nJ,1));
            intersect = false(nJ,nT);
            t = zeros(nJ,nT,'single');
            u = zeros(nJ,nT,'single');
            v = zeros(nJ,nT,'single');
            
            n0 = [0;0;1];
            q0 = [0;0;min(V(:,3))]; % horizontal plane at min(Z) to be used outside mesh
            
            for k = 1:nchunks
                
                if k*nJ <= nP
                    inj = (k-1)*nJ + 1 : k*nJ;
                else
                    inj = (k-1)*nJ + 1 : nP;
                    nJ = numel(inj);
                    
                    v1 = single(repelem(V(T(:,1),:),nJ,1));
                    v2 = single(repelem(V(T(:,2),:),nJ,1));
                    v3 = single(repelem(V(T(:,3),:),nJ,1));
                    ss(nJ*nT+1:end,:) = [];
                    intersect(nJ+1:end,:) = [];
                    t(nJ+1:end,:) = [];
                    u(nJ+1:end,:) = [];
                    v(nJ+1:end,:) = [];
                end
                vprj = single(Vp(inj,:));
                VV = repmat(vprj,nT,1);
                
                [intersect(:),t(:),u(:),v(:)] = TriangleRayIntersection(VV,ss,v1,v2,v3,options);

                ni = sum(intersect,2);
                if any(ni > 1)
                % multiple mesh intersections: pick the one closest to ray origin
                    for j = find(ni > 1)'
                        idx = find(intersect(j,:));
                        [~,tmin] = min(t(j,idx));
                        intersect(j,:) = false;
                        intersect(j,idx(tmin)) = true;
                    end
                end
                ui = u(intersect);
                vi = v(intersect);

                [pt,tr] = find(intersect);
                vprj(pt,:) = (1-ui-vi).*V(T(tr,1),:) + ui.*V(T(tr,2),:) + vi.*V(T(tr,3),:);

                outside = ni < 1;
                if any(outside)
                % get intersection with flat plane outside mesh
                    vprj(outside,:) = vprj(outside,:)-((vprj(outside,:)-q0')*n0/dot(n0,s))*s';
                end
                
                Vp(inj,:) = vprj;
            end
        end
	end
end
