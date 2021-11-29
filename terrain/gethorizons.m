function HP = gethorizons(Trck,Terr)
% HP = GETHORIZONS(TRCK,TERR) - Estimates individual near-horizon-profiles (HP.nHor) at the center
%   of each mount in TRCK, and a far-horizon-profile (HP.fHor) as the intersection of all near-
%   horizons. 
%
%   Terrain-elevation data can be provided in structure TERR, whose fields {x,y,z} represent a 
%   cloud of points over the land surface (in project coordinates! - see READTERRAIN), and
%   TERR.g0 provides the 'representative-minimum' distance among points (if most points come 
%   from a grid, e.g. a DEM raster, TERR.g0 should be the lattice of that grid). It is recommended
%   that the data provides a ~25 km margin around the site (See GETGEOTIFF for an example), 
%   with a good resolution near TRCK.centers (see READTERRAIN for merging of point clouds).
%
%   If TERR is not provided, an approximate ground mesh will be estimated by triangulating the 
%   points below the mount pivot points in TRCK (see APPROXGROUNDMESH).
%
%   In order to perform the calculations for a large number of mounts in a reasonable time, the
%   strategy used by GETHORIZONS is to divide the point cloud in several 'rings' and use a coarser
%   approximation for those which are far away from the site:
%
%   R(1) a circle that barely circumscribes the plant, and sets a lower limit for R(2)
%   R(2) near-ground. Points inside R(2) are treated as a mesh (using Delaunay-triangulation). 
%       Their contribution to near-horizon is calculated by projecting each mesh-triangle onto a
%       unit-sphere, merging the planar-stereographic projections of all triangles, and projecting
%       back to the unit-sphere.
%   R(3:n) rings beyond R(2) are designed to have at least MINRINGWIDTH points in width, 
%       (i.e. R(3) - R(2) > g0·MINRINGWIDTH) so that their horizon profile can be approximated 
%       by projecting points, not triangles. The ring-width is set so that the angular error  
%       derived from considering all points in a ring R(j-1) < d < R(j) as being projected onto 
%       a sphere R(j-1) remains below Angular Tolerance sqrt(6·RELTOL). As such, all points in 
%       each ring are projected into its 'inner-sphere', and a ring-horizon is approximated using
%       a 'moving-maximum' over a window w = 1.5·atan(g0/R(2)) - so that at least one of the
%       inner-points of the ring is considered. These ring-horizons replace all points beyond R(2)
%       and are re-projected and merged for each point of analyis.
%
% HP.nHor: size(Trck,2) set of polygon3d objects with cartesian points over the unit-sphere.
% HP.fHor: single polygon3d object with cartesian points over the unit sphere.
%
% See also: READTERRAIN, APPROXGROUNDMESH, GUITERRAIN
        
    minringwidth = getSimOption('terrain.minringwidth');
    dmin = getSimOption('terrain.minmeshdist');
    dmax = getSimOption('terrain.maxmeshdist');
    tol = getSimOption('terrain.meshtol');
    
    angtol = getSimOption('angtol')*pi/180;  % Radians! : 1 - sin(x/2)/(x/2) ~ reltol
    
    % Re-center the plant at its circumcenter
    [r0,CC] = MinBoundCircle(Trck.centers(1:2,:)');
    CC(3) = mean(Trck.centers(3,:));
    
    Trck.centers = Trck.centers - CC(:);

    if nargin > 1 && ~isempty(Terr)
    % Terrain cloud-point imported from somewhere
    
        parsestruct(Terr,{'x','y','z'},'-n','-r','-v','-e');
        parsestruct(Terr,{'g0'},'-n','-r','-s','-p');
    
        g0 = Terr.g0; % characteristic gap between points
        x = Terr.x - CC(1); 
        y = Terr.y - CC(2); 
        z = Terr.z - CC(3);
    else
    % Get an approximate mesh from the points below the mounts
        [V,~] = approxgroundmesh(Trck);
        x = V(:,1) - CC(1); y = V(:,2) - CC(2); z = V(:,3) - CC(3);
        [~,~,~,g0] = rectangleproperties(Trck.geom.border); 
        g0 = 10*g0; % ensures that all points fit in the first ring
    end
    d = hypot(x(:),y(:)); % radius of each point from site 'center'
    
    [~,~,~,L] = rectangleproperties(Trck.geom.border);
    L = L/angtol; % mount-based, minimum radius to avoid parallax error > angtol

    % Divide DEM points in 'rings':
    r = r0 + max(g0,L); % innermost circle r(1) barely circumscribes the plant
    r(2) = sqrt(g0*minringwidth/(r0*angtol)-1)*r0; % r(2), r(3) - r(2) > MIN_RING_WIDTH*g0
    if r(2) < r(1), r = r(1); end

    % get remaining rings so that ring-width doesn't produce errors > angtol
    m0 = max(d);
    while r(end) <= m0
        r(end+1) = r(end) + r(end)^2*angtol/r0; %#ok<AGROW>
    end
    
    % classify points into radius-bins, B == 1 for first circle, B == 2 for next, etc.
    r = [0,r];
    Nr = numel(r);
    [N,~,B] = histcounts(d,r);
    
    % Create a Near-Ground triangulation with points within R(2)
    idx = find(B < 3);
    T = delaunay([x(idx),y(idx)]);
    V = [x(idx),y(idx),z(idx)];
    % [V,T] = cleantrisurf(V,T,dmin,dmax,'mesh_eps1',tol/dmax,'mesh_eps2',tol/dmax);
    nGnd = triangulation(T,V);

    % Plot near-ground
    GUIfigure('horizons','Horizon - Profiles'); clf(); hold on; axh = gca();
    plothalfdome(axh,r(3),Trck.rotation);
    plotmounts(Trck);
    colormap('summer');
	trisurf(nGnd,'edgecolor','none');
    set(gca,'clipping','off');
    
    plothorizon = @(H,R,varargin) ...
        plot3(axh,H.x([1:end,1])*R,H.y([1:end,1])*R,H.z([1:end,1])*R,varargin{:});
 
    if Nr > 3
    % Calculate mid-range (ring) horizons
        set(gca,'colororder',summer(Nr-3));
        waitbarh = optwaitbar(0,'Getting ring-horizons','Name','Calculating Far-Horizon Profiles...');

        w = 1.5*atan(g0/r(3)); % moving window width
        mHor(Nr-3) = polygon3d();
        for j = 3:Nr-1
            idx = find(B == j);
            [az,el] = cart2sph(x(idx),y(idx),z(idx));
            %[az,ic] = sort(az); el = el(ic); 

            % DEBUG
            %if ishandle(h), delete(h); end
            % Vp = [x(idx),y(idx),z(idx)];
            %Vp = bsxfun(@rdivide,Vp,rssq(Vp,2));
            %Vp = stereoproj(Vp(:,1),Vp(:,2),Vp(:,3));
            %scatter(Vp(:,1),Vp(:,2),2);
            %scatter3(Vp(:,1),Vp(:,2),Vp(:,3),5);

            th = 0:angtol:2*pi;
            ph = zeros(size(th));
            src = zeros(size(th)); % point source
            for k = 1:numel(th)
                angdiff = rem(az - th(k) + 3*pi, 2*pi) - pi; % angle difference, limited to ±180°
                idx = abs(angdiff) < w/2;
                if ~any(idx)
                    src(k) = -k; % just mark it as unique*
                    continue; 
                end
                [ph(k),ic] = max(el(idx));
                idx = find(idx);
                src(k) = idx(ic);
                th(k) = az(src(k));
                waitbarh.update((sum(N(3:j-1))+ N(j)*k/numel(th))/sum(N(3:end)),...
                    sprintf('Ring %d / %d : %0.1f to %0.1f km',j,Nr-1,r(j)/1000,r(j+1)/1000));
            end
            [~,ic] = unique(src,'stable');
            th = th(ic);
            ph = ph(ic);
            
            % Project onto a shpere with min. ring radius (nearest objects have less offset)
            mHor(j-2) = polygon3d(th*180/pi,ph*180/pi,r(j),'sph');
            
            % Plot just outside r(3)- sphere
            plothorizon(mHor(j-2),r(3)*(0.8+j/10)/r(j));
            drawnow();
        end
        delete(waitbarh);
    else
        mHor = polygon3d.empty;
    end
    
    % Get actual points of analysis (considering rotated non-zero axisoffset)
    % Currently only centroids of every tracker, at a single solar position...
    % FUTURE: calculate horizon on-the-fly, use multiple viewpoints
    pts = getanalysedpts(Trck);

    % Calculate local-horizons
    M = getnworkers([],'vars',whos('nGnd','mHor','prj')); 
    if M > 1
        HP = runparallel(@localhorizons,{pts,nGnd,mHor,angtol},1);
    else
        HP = localhorizons(pts,nGnd,mHor,angtol);
    end
%     if iscell(HP.fHorP)
%         HP.fHorP = mergepolygons(polygon.empty,[HP.fHorP{:}],'pos','pos');
%     else
%         HP.fHorP = mergepolygons(polygon.empty,HP.fHorP,'pos','pos');
%     end          
%     HP.fHorP([HP.fHorP.hole]) = [];
%     HP.fHor = prj.inverse(HP.fHorP);
    
    % Plot near horizon
    c = cool(size(pts,1));
    arrayfun(@(h,j) plothorizon(h,r(3),'color',c(j,:)),HP.nHor,(1:size(pts,1))');
    plothorizon(HP.fHor,r(3),'k-','linewidth',2);
end

function pts = getanalysedpts(Trck)
% Currently only centroids of every tracker, at a single solar position...
% FUTURE: calculate horizon on-the-fly, use multiple viewpoints

    c = centroid(Trck.geom.border);
    Trck.geom = pvArea();
    [Trck.geom.border.x,Trck.geom.border.y] = deal(c(1),c(2));
    [~,pts] = plantlayout(Trck,0,90-abs(Trck.origin(2)),'Eq0');
end

function varargout = plotmounts(Trck)
% Minimal Layout-plot (only mount patches, sun at equinox noon)

    [F,V] = plantlayout(Trck,0,90-Trck.origin(2),'eq0');
    h = patch('Faces',F,'Vertices',V,'FaceColor',[0.3,0.3,0.8],'edgecolor','w','edgealpha',0.1);
    
    if nargout > 0, varargout{1} = h; end
end

function HP = localhorizons(pts,nGnd,mHor,angtol)
% Get horizons for points PTS, given triangulation nGnd, and horizon profiles
        
    % validateattributes(pts,{'numeric'},{'real','size',[NaN,3]},'','PTS');
    % validateattributes(nGnd,{'triangulation'},{'nonempty'},'','nGnd');
    % validateattributes(nGnd.Points,{'numeric'},{'real','size',[NaN,3]},'','nGnd.Points');

    waitbarh = optwaitbar(0,'Getting ring-horizons','Name','Calculating Near-Horizon Profiles...');
    np = size(pts,1);
        
    T = nGnd.ConnectivityList;
    
    % Create a stereographic projector, to merge near-triangles on a plane
    prj = polyprojector('stereo','clip',179.9);
    
    % Flip everything upside down, so that the nadir maps to 0,0
    V = nGnd.Points;
    V(:,3) = -V(:,3);
    pts(:,3) = -pts(:,3);
    for k = 1:numel(mHor)
       mHor(k).z = -mHor(k).z;
    end
    
    % Assume the horizon does not cover the zenith*
    % (* cannot use [0,0,1]' because it maps to infinity)
    OUTSIDE_PT = [1e-3;0;-1];
    
    nHor(np,1) = polygon3d();
    areaerror = zeros(size(pts,1)+1,1);
    
    fHorP = polygon(1:360,2,'pol');
    
    for j = 1:size(pts,1)
        pt = pts(j,:)';

        % start with a flat horizon
        L = polygon(1:360,1,'pol');
        if ~isempty(mHor)
        % Project ring-horizons
            q = arrayfun(@(p) prj.project(polytranslate(p,-pt),OUTSIDE_PT,false),mHor);
            L = mergepolygons(L,q,'pos','pos');
        end

        % Stereographic projection of all vertices
        Vo = V-pt';
        Vp = prj.prj(Vo');

        % Note which triangles lie completely below the horizon
        vbh = insidepolygon(L,Vp(1,:),Vp(2,:));   % vertices below horizon
        tbh = all(vbh(T),2);                      % triangles below horizon

        if ~all(tbh)
        % Merge all (~tbh) triangle projections with the current-horizon
            P = polygon3d.vf2poly(Vo,T(~tbh,:));
            P = prj.project(P,OUTSIDE_PT,false);
            % P = polygon3d.vf2poly(Vp',T(~tbh,:));
            % P = fixorientation(P);
            L = mergepolygons(L,P,'pos','pos');       % merge projected obstacles
        end
        L(abs([L.signedarea]) < (angtol^2)/6) = []; % remove splinter polygons

        if ~isscalar(L)
            % save('gethorizons_dump.mat');  % DEBUG
            a = [L.signedarea];
            [~,biggest] = max(a);
            L = L(biggest);
            areaerror(j) = sum(abs(a(~biggest)))/pi;
        end
        
        % Keep track of far-horizon as the intersection of all near-horizons
        fHorP = intersectpolygons(fHorP,L); % i.e. sky that is visible from at least one pt

        % Get back x,y,z coordinates on unit sphere from Stereographic projection
        % nHor(j) = polygon3d(invstrproj(P.x,P.y));
        nHor(j) = prj.inverse(L);
        nHor(j).z = -nHor(j).z;

        waitbarh.update(j/np,sprintf('Position %d / %d',j,np),'-addtime');

        % Remove unnecessary vertices
        % nHor(j) = cleanpolygon(nHor(j),1-cosd(angtol));
    end
    
    if ~isscalar(fHorP)
        % save('gethorizons_dump.mat');  % DEBUG
        a = [fHorP.signedarea];
        [~,biggest] = max(a);
        fHorP = fHorP(biggest);
        areaerror(end) = sum(abs(a(~biggest)))/pi;
    end
    
    fHor = prj.inverse(fHorP);
    fHor.z = -fHor.z;
    
    if max(areaerror) > 0
        warning(['Failed to find %d continuous horizon profile(s). Dumping splinter polygons '...
                 'with up to %0.2f %% of dome)'],nnz(areaerror),max(areaerror));
    end
    
    HP = struct('nHor',nHor,'fHor',fHor);
end