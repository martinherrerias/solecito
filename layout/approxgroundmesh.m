function varargout = approxgroundmesh(Mounts,buffer,dmin,dmax,tol)
% [V,T] = APPROXGROUNDMESH(MOUNTS,...)
% TRI = APPROXGROUNDMESH(MOUNTS,[BUFFER,DMIN,DMAX]) - where TRI = triangulation(T,V)
%
%   Generate a TIN mesh guessing what the ground might look like beneath a set of MOUNTS, 
%   based on the set of coordinates and dimensions of MOUNTS. For 2a, 1aV, and 1aF mount types, 
%   (see MOUNTROTATIONS for definitions) this involves using a single point per mount, a 
%   distance MOUNTS.centerheight below the provided MOUNTS.centers. For 0a and 1aC, it is 
%   assumed that the height of the terrain at the x-axis edges (!), and not on the middle point,
%   determines the mount position. The mounts are therefore rotated (to solar noon), their edge
%   positions calculated, and these points (minus MOUNTS.centerheight on Z) are used in the
%   triangulation.
%
%   Unless the BUFFER flag is set to false (default is true), the triangulation will also include
%   a buffer polygon with a (~one-mount-diagonal) margin around the convex hull of all mounts, 
%   whith elevations based on the nearest neighboring point.
% 
%   DMIN,DMAX are used to simplify/refine the mesh obtained as described above. Defaults are
%   taken from SimOptions.terrain.minmeshdist/maxmeshdist.
%
% [TRI,VV,P,w] = APPROXGROUNDMESH(MOUNTS,...) call FINEMESHPOINTS with n = floor(DMAX/DMIN) to
%   provide a finer set of vertices VV, of which P(k,:) are located in TRI(k,:).
%
% See also: PLOTTRACKERARRAY, READTERRAIN

    SHRINK_FACTOR = 0.5;

    parsestruct(Mounts,{'centers'},'numeric','real','finite','nonempty','size',[3,NaN]);
    parsestruct(Mounts,{'centerheight'},'numeric','real','scalar','positive');
    parsestruct(Mounts,{'axisoffset'},'numeric','real','finite','vector','numel',3);
    parsestruct(Mounts,{'geom'},'class','pvArea','scalar',@(x) isscalar(x.border));
    Mounts.type = parselist(Mounts.type,{'2a','1aV','1aF','1aC','0a'});
    
    [w,h] = rectangleproperties(Mounts.geom.border); % tracker width/height

    if nargin < 2, buffer = 0; end
    
    if nargin < 3 || isempty(dmin), dmin = getSimOption('terrain.minmeshdist'); end
    if nargin < 4 || isempty(dmax), dmax = getSimOption('terrain.maxmeshdist'); end
    if nargin < 4 || isempty(tol), tol = getSimOption('terrain.meshtol'); end
    validateattributes(dmin,{'numeric'},{'scalar','positive','real','finite'});
    validateattributes(dmax,{'numeric'},{'scalar','positive','real','>=',dmin});
    dmin = min([w/2,h/2,dmin]);
    dmax = max([dmax,2*dmin]);
    
    V = getbasepoints(Mounts,dmin);
    
    if isequal(buffer,'raw'), varargout{1} = V; return; end
            
    validateattributes(buffer,{'numeric'},{'scalar','nonnegative','real','finite'});
    buffer = max([buffer,dmax]);

    [V,T] = boundedsurfacetriangulation(V,buffer,SHRINK_FACTOR);
   
    % V0 = V; T0 = T;
    [V,T] = cleantrisurf(V,T,dmin,dmax,'mesh_eps1',tol/dmax,'mesh_eps2',tol/dmax);
    
    % GUIfigure('debug'); clf();
    % subplot(1,2,1);
    % triplot(T,V(:,1),V(:,2)); hold on;
    % axis equal; %  axis([-5 20 10 40])
    % 
    % TR = triangulation(T,V(:,1:2));
    % [ID,B] = TR.pointLocation(V0(:,1:2));
    % outside = isnan(ID);
    % plot(V0(outside,1),V0(outside,2),'ro');
    % Z = V(:,3);
    % z = dot(Z(TR(ID(~outside),:)),B(~outside,:),2);
    % subplot(1,2,2);
    % hist(V0(~outside,3)-z);
  
    if nargout == 2, [varargout{1:2}] = deal(V,T);
    else
        T = triangulation(T,V);
        varargout{1} = T;
        if nargout > 1
            E = T.edges;
            dmax = prctile(rssq(V(E(:,1),:)-V(E(:,2),:),2),75);
            n = floor(dmax/dmin);
            [varargout{2:4}] = finemeshpoints(varargout{1},n);
            
        end
    end
end

function V = getbasepoints(Mounts,dmin)
    % get tracker base points
    x = Mounts.centers(1,:);
    y = Mounts.centers(2,:);
    z = Mounts.centers(3,:)-Mounts.centerheight;

    % For 1aC and 0a mounts, use axis pivot points, instead of centers
    if contains(Mounts.type,{'1aC','0a'})

        RotMat = mountrotations(Mounts,0,90,'N2E');
        Ntr = size(Mounts.centers,2);

        % Get axis end-point vectors
        ax0(:,1) = [max(Mounts.geom.border.x);0;0] + Mounts.axisoffset(:);    
        ax0(:,2) = [min(Mounts.geom.border.x);0;0] + Mounts.axisoffset(:);

        % Rotate them...
        if size(RotMat,3)==1
            offsets(:,1,1) = RotMat(:,:,1,1)*ax0(:,1);
            offsets(:,1,2) = RotMat(:,:,1,1)*ax0(:,2);
            offsets = repmat(offsets,[1,Ntr,1]);
        else
            for j = Ntr:-1:1
                offsets(:,j,2) = RotMat(:,:,j,1)*ax0(:,2);
                offsets(:,j,1) = RotMat(:,:,j,1)*ax0(:,1);
            end
        end

        % Replace centers with pivot points
        x = [x+offsets(1,:,1),x+offsets(1,:,2)];
        y = [y+offsets(2,:,1),y+offsets(2,:,2)];
        z = [z+offsets(3,:,1),z+offsets(3,:,2)];
    end
    
    [V,~,ia] = uniquetol([x;y]',dmin/2,'DataScale',1,'ByRows',true);
    V(:,1) = accumarray(ia,x',[],@mean);
    V(:,2) = accumarray(ia,y',[],@mean);
    V(:,3) = accumarray(ia,z',[],@mean);

%     V = [x;y;z]';

%     % Remove redundant points
%     TooClose = pdist(V) < dmin;
%     % TC = squareform(pdist([x',y',z']) < DMIN);
%     
%     M = numel(x);
%     tooclose = false(1,M); 
%     keepers = true(1,M);
%     for j = 1:M
%         if ~keepers(j), continue; end
%         tooclose(:) = false;
%         tooclose(j+1:M) = TooClose((j-1)*(M-j/2)+(1:M-j)) & keepers(j+1:M);
%         % tooclose = keepers & TC(j,:);
%         if any(tooclose)
%             tooclose(j) = true;
%             V(j,:) = mean(V(tooclose,:),1);
%             keepers(tooclose) = false;
%             keepers(j) = true;
%         end
%     end
%     V = V(keepers,:);
end

function [V,T] = boundedsurfacetriangulation(V,buffer,SHRINK_FACTOR)

    if rank(V) < 2
        p = polygon.outline(V(:,1),V(:,2),buffer,6);
    else
        k = boundary(V(:,1),V(:,2),SHRINK_FACTOR);
        % k = convhull(V(:,1),V(:,2));
        p = polygon(V(k,1),V(k,2));
        p = offsetpolygon(p,buffer);
    end
    border = poly2vef(refinepoly(p,buffer));
    
    [S,L] = bounds(V(:,1:2),1);
    r = hypot(diff(S),diff(L)); 
    r = max(r,buffer)*1000; % very far away
    exterior = poly2vef(offsetpolygon(p,r)); 
    
    Z = V(:,3);
    
    [w,iz] = pdist2(V(:,1:2),border,'euclidean','Smallest',100);
    w = max(eps,w').^(-4);
    w = w./sum(w,2);
    border(:,3) = dot(Z(iz'),w,2);

%     % Use inverse-square interpolation to get heights far away
%     [~,iz] = pdist2(V(:,1:2),exterior,'euclidean','Smallest',2);
%     Z = [Z;Z(iz)];
    
%     DT = delaunayTriangulation([V(:,1:2);exterior]);
%     [ID,B] = DT.pointLocation(border);
%     border(:,3) = dot(Z(DT(ID,:)),B,2);
    
    V = [V;border];
    
    T = delaunay([V(:,1:2);exterior]);
    T(any(T > size(V,1),2),:) = [];
end
