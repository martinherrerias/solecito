function [V,K,B,s] = voronoisphere(xyz,varargin)
% [V,K,B,s] = VORONOISPHERE(XYZ) - Compute the Voronoi diagram of points XYZ on the S² sphere 
%
% INPUT:
%   XYZ - (3 x n) array, coordinates of n points in R^3
%       NOTE: they are expected to be normalized to 1 (i.e., belong to the 2-sphere) and unique.
%
% OUTPUT:
%   V - (3 x m) array, coordinates of the vertices of voronoi diagram
%   K - (n x 1) cell, each K{j} has indices of the voronoi-cell vertices for xyz(:,j)
%       Vertices are counter-clockwise oriented when looking from outside.
%   B - Boundary, (n x 1) cell, each B{j} contains the discretized spherical polygonal coordinates
%       of the voronoi cell. They are (3 x mj) arrays, mj are number of discretized points.
%   s: (n x 1) array, solid angle of voronoi cell
%
% VORONOISPHERE(..., 'resolution', rrad) to provide the resolution of the
%   boundary. RRAD is in radian. Default value is 0.0349 (~ 2 degrees).
%   
% Original Author: Bruno Luong <brunoluong@yahoo.com> 
% Date creation: 28/March/2013
%                29/March/2013: specify orientation
%                17/July/2018: handle small number of data (<=3), output solid angle
%
% Modified: Martín Herrerías <martin.herrerias@gmail.com>
%           23/Aug/2018: improve documentation
%                        simplify using MATLAB's TRIANGULATION class
%                        fix issues with flat hulls and multiple points in great-arc.
%
% See also: VORONOI, VORONOIN, TRIANGULATION

%{
Original Copyright (c) 2013, Bruno Luong
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution
* Neither the name of FOGALE nanotech nor the names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%}

%{
% Example
n = 100;
xyz = randn(3,n);
xyz = bsxfun(@rdivide, xyz, sqrt(sum(xyz.^2,1)));

[P, K, B] = voronoisphere(xyz);

%% Graphic
f = figure(1);
clf(f);
set(f,'Renderer','zbuffer');
ax = axes('Parent', f);
hold(ax, 'on');
set(ax, 'Color', 'w');
plot3(ax, xyz(1,:),xyz(2,:),xyz(3,:),'wo');
clmap = cool();
ncl = size(clmap,1);
for k = 1:n
    X = B{k};
    cl = clmap(mod(k,ncl)+1,:);
    fill3(X(1,:),X(2,:),X(3,:),cl,'Parent',ax,'EdgeColor','w');
end
axis(ax,'equal');
axis(ax,[-1 1 -1 1 -1 1]);
%}

    MAXSIZE = 2500;

    if size(xyz,1) ~= 3, xyz = xyz'; end
    assert(size(xyz,1) == 3,'Expecting 3xN array');
    
    xyz = xyz./rssq(xyz,1);
    if size(xyz,2) > MAXSIZE
        assert(isreasonable(size(xyz,2)),'Too many points!');
    end

    % Get the resolution in radianS
    [opt] = getpairedoptions(varargin,{'resolution'},{2*pi/180}); % 2 degrees
    resolution = opt.resolution;
    
    casetype = size(xyz,2);
    if rank(xyz - mean(xyz,2)) < 3, casetype = 3; end   % all points on a single plane
    
    switch casetype
    case 0
        V = zeros(3,0);
        K = cell(0,1);
        B = cell(0,1);
        if nargout >= 4
            s = zeros(0,1);
        end
    case 1
        V = zeros(3,0);
        K = {zeros(0,1)};
        B = {zeros(3,0)};
        if nargout >= 4
            s = 4*pi;
        end      
    case 2
        V = zeros(3,0);
        K = repmat({zeros(0,1)},[2 1]);
        B{1} = FullCircle(diff(xyz,1,2), resolution);
        B{2} = fliplr(B{1});
        if nargout >= 4
            s = 2*pi+zeros(2,1);
        end
    case 3
        % If all points lie on a single plane...
        c = mean(xyz,2);
        Q = orth(xyz-c);                     % get reduced orthonormal basis (3x2 matrix)
        if size(Q,2) > 2, Q = Q(:,1:2); end  % ... (fix precision issues)
        xy = Q'*(xyz-c);                	 % ... project to plane
        idx = convhull(xy');
        xyz = xyz(:,idx);                    % ... and reorder points around flat convex-hull
   
        nf = size(xyz,2)-1;                  % (convhull repeats first point at the end)
        dxyz = diff(xyz,1,2);
        
        % All triangles share a common center (normal to plane)
        V = cross(Q(:,2),Q(:,1));

        % Build the great-arcs that link [-V,V] perpendicular to edges
        edgearcs = arrayfun(@(k) HalfCircle(V, dxyz(:,k), resolution),1:nf,'unif',0);

        % Build the contour of the voronoi cells
        B = cell(nf,1);
        for k = 1:nf
            B{k} = cat(2, edgearcs{mod(k-2,nf)+1}, edgearcs{k}(:,end-1:-1:2));
        end
        
        % All cells just go (strictly speaking) from one pole to the other
        K = repmat({[1 2]},[1 nf]);
        
        % Revert to original order
        [~,idx] = sort(idx(1:end-1));    
        K = K(idx);
        B = B(idx);

        if nargout >= 4
            dxyz_p1 = dxyz(:,[2:end,1]);
            a = dxyz./rssq(dxyz,1);
            b = bsxfun(@cross,a,V);
            s = 2*atan2(dot(dxyz_p1,b),dot(dxyz_p1,a));
            
            s = s(idx);
        end

        % % DEBUG: Plot triangulation, straight edges and arcs
        % clf; hold on
        % arrayfun(@(j) plot3(xyz(1,[1:end,1]),xyz(2,[1:end,1]),xyz(3,[1:end,1]),'r:'),1:nf+1);
        % arrayfun(@(j) text(xyz(1,j),xyz(2,j),xyz(3,j),num2str(j)),1:size(xyz,2));
        % arrayfun(@(j) text(V(1,j),V(2,j),V(3,j),num2str(j),'color','g'),1:size(V,2));
        % cellfun(@(b) plot3(b(1,:),b(2,:),b(3,:),'c-'),B);
        % axis equal;

    otherwise
        T = convhull(xyz.');
        TR = triangulation(T,xyz');
        tri_edges = TR.edges;

        % For each triangulation-edge j, there is a perpendicular Voronoi-cell-edge that 
        %   joins the two circumcenters [p q] = cell_edges(j,:)
        cell_edges = cell2mat(TR.edgeAttachments(tri_edges));

        % Center of the circumscribed Delaunay triangles
        V = circumcenters(xyz,T);

        %TOL = min(pdist(xyz'))/10000;
        %TOL = max([1e-12,1/size(xyz,2),TOL]);
        TOL = 1E-8;
        
        % Check that triangle centers do not match (e.g. two triangles on same great-arc)
        [~,ia,ic] = unique(round(V'/TOL),'rows');
        if numel(ia) ~= numel(ic)
            V = V(:,ia);
            cell_edges = ic(cell_edges);
            degenerate = cell_edges(:,1) == cell_edges(:,2);
            cell_edges = cell_edges(~degenerate,:);
            tri_edges = tri_edges(~degenerate,:);
        end

        % Collect all the edges that share a given triangulation vertex into each Voronoi-cell:
        %   faces{j} contains a list of the edgeIDs (i.e. the row number in tri_edges or
        %   cell-edges) that ar
        ne = size(tri_edges,1);
        edgeIDs = repmat((1:ne).',[2 1]);
        faces = accumarray(tri_edges(:),edgeIDs, [], @(e) {e});
        nf = numel(faces);

        % Build the great-arcs that link every pair of circumcenters
        cell_edges = sort(cell_edges,2);
        edgearcs = arrayfun(@(j) Arc(V(:,cell_edges(j,:)),resolution),1:ne,'unif',0);

        % Build the contour of the voronoi cells          
        B = cell(nf,1);
        K = cell(nf,1);
        for k = 1:nf
        % ordering and orientation of the edges
            fce = cycling_edges(faces{k}, cell_edges);
            if ~oriented_ccw(V(:,fce([1:end,1],1)),xyz(:,k))
                fce = rot90(fce,2);
            end
            K{k} = fce(:,1); % Keep the indices of the voronoi's hull

            % joint the corresponding arcs to form the Boundary
            [~, loc] = ismember(sort(fce,2), cell_edges, 'rows');
            X = edgearcs(loc);
            flip = fce(:,1) ~= cell_edges(loc,1);
            X(flip) = cellfun(@fliplr, X(flip), 'unif', false);
            X = cellfun(@(x) x(:,1:end-1), X, 'unif', false);   % remove duplicat end-points
            B{k} = cat(2, X{:});
        end

        if nargout >= 4
            s = vcell_solidangle(V, K, xyz);
        end

        % % DEBUG: Plot triangulation, straight edges and arcs
        % clf; hold on
        % arrayfun(@(j) plot3(xyz(1,tri_edges(j,:)),xyz(2,tri_edges(j,:)),xyz(3,tri_edges(j,:)),'r:'),1:size(tri_edges,1));
        % arrayfun(@(j) plot3(V(1,cell_edges(j,:)),V(2,cell_edges(j,:)),V(3,cell_edges(j,:)),'y-'),1:size(cell_edges,1));
        % arrayfun(@(j) text(xyz(1,j),xyz(2,j),xyz(3,j),num2str(j)),1:size(xyz,2));
        % arrayfun(@(j) text(V(1,j),V(2,j),V(3,j),num2str(j),'color','g'),1:size(V,2));
        % cellfun(@(b) plot3(b(1,:),b(2,:),b(3,:),'c-'),B);
        % axis equal;
    end
end % voronoisphere

function G = FullCircle(Pole, resolution)
% Return an full circle perpendicular to Pole
    Pole = Pole/norm(Pole);
    Q = null(Pole.');
    A = Q(:,1);
    B = cross(A,Pole); % Q(:,2);
    theta = 2*pi;
    npnts = max(ceil(theta/resolution),2);
    theta = linspace(0, theta, npnts);
    G = A*cos(theta) + B*sin(theta);
end

function G = HalfCircle(A, Pole, resolution)
% return in G a half circle that links A to -A perpendicular to Pole
    Pole = Pole/norm(Pole);
    A = A/norm(A);
    B = cross(A,Pole);
    theta = pi;
    if nargin >= 3
        npnts = max(ceil(theta/resolution),2);
        theta = linspace(0, theta, npnts);
        G = A*cos(theta) + B*sin(theta);
    else
        G = zeros(3,0);
    end
end

function G = Arc(AB, resolution)
% Return a discretized arc between points A and B
    A = AB(:,1);
    B = AB(:,2);
    AxB = cross(A,B);
    AdB = dot(A,B);
    Ap = cross(AxB, A);
    Ap = Ap/norm(Ap);
    theta = atan2(sqrt(sum(AxB.^2,1)), AdB); % > 0
    npnts = max(ceil(theta/resolution),2); % at least 2 points
    theta = linspace(0, theta, npnts);
    G = A*cos(theta) + Ap*sin(theta);
end % Arc

function P = circumcenters(xyz,T)
% Centers of the circumscribed Delaunay triangles for triangulation T on vertices XYZ
%
% In most cases, the fuction is (within eps) equivalent to:
%
%   TR = triangulation(T,xyz');
%   C = TR.circumcenter';
%   P = C./rssq(C,1);
%
% However, when all points in a given triangle fall over a great-arc, C gets to [0,0,0] and this
% no longer works.

    XYZ = reshape(xyz(:,T),[3 size(T)]);
    A = XYZ(:,:,1);
    B = XYZ(:,:,2);
    C = XYZ(:,:,3);
    A = A-C;
    B = B-C;
    A2B = bsxfun(@times, sum(A.^2,1), B);
    B2A = bsxfun(@times, sum(B.^2,1), A);
    AxB = cross(A,B,1);
    P = cross(A2B - B2A, AxB, 1);
    P = C + bsxfun(@times,P,1./(2*sum(AxB.^2,1)));
    nP = sqrt(sum(P.^2,1));
    P = bsxfun(@times, P, 1./nP);
    s = dot(AxB,C);
    P = bsxfun(@times, P, sign(s));

    if any(s==0)
        P(:,s==0) = AxB(:,s==0)./rssq(AxB(:,s==0));
    end
end

function C = cycling_edges(F,E)
% Given a set of edges F (1xN array of edge-indices) and the definition E of those edges (Mx2 array 
% of vertex-indices, chain the edges in an Nx2 cycle of vertices C = [a,b;b,c;c,d;...;x,a]

    u = E(F,:).';
    n = size(u,2);
    [~,~,I] = unique(u);
    I = reshape(I,[2 n]);  % index of vertices, edge j contains I(:,j)
    J = repmat(1:n,[2 1]);
    if ~all(accumarray(I(:), 1) == 2)
        error('Topology issue due to numerical precision')
    end
    K = accumarray(I(:), J(:), [], @(x) {x}); % to which edges K(:,j) does vertex j belong?
    K = [K{:}];
    
    % chain the edges
    C = zeros([n 2]);   % sequence of vertices
    q = 1;              % current vertex
    p = K(2);           % current edge - initialize at K(2,1) makes next-edge K(1,1) (*)
    for j = 1:n
        i = K(:,q);     % there are 2 edges connected to vertex q ...
        p = i(i ~= p);  % next edge is the one which is not p                        (*)
        i = I(:,p);     % there are 2 vertices in new edge p ...
        if i(1) == q
        % if vertex order is already ok, move to second vertex
            C(j,:) = u([1 2],p);
            q = i(2);
        else
        % otherwise flip vertices, then move to next
            C(j,:) = u([2 1],p);
            q = i(1);
        end
    end
end

function ccw = oriented_ccw(V,c)
% Make sure the set of vertices V [3xN] is oriented counter-clockwise around point c [3x1]

    Q = null(c.');              % orthogonal components to c
    xy = Q'*V;                  % flat vertex projection (perpendicular to c)
    a = (xy(1,1:end-1)-xy(1,2:end))*(xy(2,1:end-1)+xy(2,2:end))'; % signed area
    
    % Combine orientation and 'directness' of Q to check for CCW orientation
    % (the sign of the determinant determines if a linear transformation preserves orientation)
    ccw = xor(a < 0, det([c Q]) < 0);
    
end

function s = vcell_solidangle(P, K, xyz)
% s = vcell_solidangle(P, K)
% s = vcell_solidangle(P, K, xyz)
%
% Compute the solid angles of All voronoi cell
%
% P is (3 x m) array of vertices, coordinates of the vertices of voronoi diagram
% K is (n x 1) cell, each K{j} contains the indices of the voronoi cell
% xyz is (3 x n) optional knot points to guide vcell_solidangle to compute
%   the solid angle of the "right" cell containing the node (and not the
%   complement cell)
%
% Restrictions:
% - P must be unit vectors
% - For each cell vertices must be counter-clockwise oriented when looking from outside
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% Date creation: 16/July/2014
%
% See also: voronoisphere

    % Turn it to false to maximize speed
    check_unit = true;

    if check_unit
        u = sum(P.^2, 1);
        if any(abs(u-1) > 1e-6)
            error('vcell_solidangle:Pnotnormalized', ...
                'vcell_solidangle: P must be unit vectors');
        end
    end

    if nargin < 3
        s = cellfun(@(k) one_vcell_solidangle(P(:,k)), K);
    else
        n = size(xyz,2);
        s = arrayfun(@(j) one_vcell_solidangle(P(:,K{j}),xyz(:,j)), (1:n).');
    end
end % vcell_solidangle

function omega = one_vcell_solidangle(v, center)
% omega = one_vcell_solidangle(v)
% Compute the solid angle of spherical polygonal defined by v
% v is (3 x n) matrix, each column is the coordinates of the vertex (unit vector, not checked)
% "correctly" oriented
%
% Ref: A. van Oosterom, J. Strackee: "A solid angle of a plane triangle."
%   IEEE Trans. Biomed. Eng. 30:2 (1983); 125-126. 
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% Date creation: 16/July/2014
%
% See also: vcell_solidangle

    if nargin < 2
        n = size(v,2);
        s = zeros(1,n-2);
        for i = 2:n-1
            T = v(:,[1 i i+1]);
            num = det(T);
            denom = 1 + sum(sum(T .* T(:,[2 3 1]), 1), 2);
            s(i-1) = num / denom;
        end
    else
        v(:,end+1) = v(:,1);
        n = size(v,2);
        s = zeros(1,n-1);
        for i = 1:n-1
            T = [center, v(:,[i i+1])];
            num = det(T);
            denom = 1 + sum(sum(T .* T(:,[2 3 1]), 1), 2);
            s(i) = num / denom;
        end
    end

    omega = atan(s);
    omega = 2 * sum(omega);

end % one_vcell_solidangle

function isit = isreasonable(N)
% Estimate the time that would take to calculate VORONOISHPERE(X) for size(X,2) = N, and ask the
% user if it seems reasonable.

    hms = @(x) sprintf('%02d:%02d:%02d',floor([x,rem(x,[3600 60])]./[3600 60 1]));

    n = round(logspace(1,min(3,log10(N)-1),10));
    for j = numel(n):-1:1
       tic();
       voronoisphere(randn(n(j),3));
       t(j) = toc();
    end
    [p,S] = polyfit(n,t,2);
    [T,dt] = polyval(p,N,S);
    msg = sprintf(['Evaluating VORONOISHPERE for %d points might take %s ± %s (HH:MM:SS), ',...
        'are you sure you want to wait?... it might just crash'], N,hms(T),hms(2*dt));
    
    switch questdlg(msg,'voronoishpere','Quit','Wait','Quit')
        case 'Wait'
            isit = true;
        otherwise
            isit = false;
    end
end