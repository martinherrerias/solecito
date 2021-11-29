function Terr = readterrain(Trck,filenames)
% TERR = READTERRAIN(TRCK,FILENAMES,OPT)
% Read a (set of) Digital-Elevation-Model file(s), and convert to a merged point-cloud TERR in 
% Project-coordinates, verifying that point coordinates fit those of the Project (TRCK), and 
% smoothing-out small differences between elevation data from the file(s) and from TRCK.centers. 
% In order for things to work smoothly:
%
%   TRCK: must provide the means of converting mount-coordinates (TRCK.centers) to EPSG 4326
%     (see CHECKCOORDSYSTEM). Additional fields are required by APPROXGROUNDMESH.
%
%  FILENAMES: READTERRAIN currently takes up to one GEOTIFF (*.tif) and one POINT-CLOUD (*.xyz)
%     files, passed as a cell-array of string file-names (if omitted, *.tif and *.xyz files will 
%     be searched-for in the working diretory).
%
%     GEOTIFF file is expected to be a raster of elevation points (meters) in WGS84 EPSG 4326 
%       (Lat,Long) coordinates. It's recommended that it provides a ~25 km margin around the 
%       project, with moderate (~1") resolution.
%
%     POINT-CLOUD file is expected to be a plain-text file with three columns: X;Y;Z.
%       X,Y coordinates should be WGS84 EPSG 4326 decimal degrees with suficient precision.
%       Z should be in meters. An optional header should be commented #X;Y;Z.
%       Points should cover the near-ground (~50 m around the project) with ~5 m resolution.
%
% Output:
%   TERR: structure with fields {x,y,z} representing both the points in the DEM(s) and the ground
%     just below TRCK (see APPROXGROUNDMESH), for better resolution close to the mounts.
%     TERR.g0 contains the mean lattice of the DEM, for use in GETHORIZONS.
%
% See also: GETHORIZONS, GUITERRAIN, APPROXGROUNDMESH

    if nargin < 2 || isempty(filenames), filenames = pickfile({'*.tif';'*.xyz'},Inf); end
    assert(~isempty(filenames),'readterrain:nofile','Need at least one file to work');
    if ~iscell(filenames), filenames = {filenames}; end
    assert(iscellstr(filenames),'Expecting string or cell-array of string filename(s)'); %#ok<ISCLSTR>
    
    [R0,C] = MinBoundCircle(Trck.centers(1:2,:)'); % characteristic (horizontal) length of plant
    
    opt.angtol = getSimOption('angtol');
    opt.buffer = getSimOption('terrain.buffer');
    opt.rmin = getSimOption('terrain.demradius');
        
    % min. DEM radius for meaningful near-horizons: anything above can be considered 'far-horizon'
    % i.e. it is the same (within angtol) for all points on the plant (see GETHORIZONS).
    opt.rmin = max(opt.rmin,R0/tand(opt.angtol)); 
    
    [~,~,~,SYS] = checkcoordsystem(0,0,0,Trck,'output','prj');
    
    % Start with pts representing mount-positions (and their inmediate ground)
    [pts,~] = approxgroundmesh(Trck);
    
    % If there is a point-cloud file available, attach the points to the list
    xyzfile = filenames(cellfun(@(x) strcmpi(x(end-2:end),'xyz'),filenames));
    if ~isempty(xyzfile)
        fprintf('Reading point-cloud...\n');
        pts = readpointcloud(xyzfile{1},pts,SYS);
    end
    g0 = mediandist(pts(:,1),pts(:,2));

    % Finally, if there is a GEOTIFF file, add it to the list of points
    gisfile = filenames(cellfun(@(x) strcmpi(x(end-2:end),'tif'),filenames));
    if ~isempty(gisfile)
        fprintf('Reading geotiff file...\n');
        [pts,g0] = readgeotiff(gisfile{1},pts,SYS,opt);
    end
    
    hfig = GUIfigure('terrain','Terrain'); clf(hfig);
    near = rssq(pts(:,1:2)-C,2) < 5*(R0 + opt.buffer);
    [~,DT] = polygon.outline(pts(near,1),pts(near,2));
    trisurf(triangulation(DT.ConnectivityList,[DT.Points,pts(near,3)]),'edgealpha',0.1);
    hold on; axis equal;
    colorbar();
    % scatter3(Trck.centers(1,:)',Trck.centers(2,:)',Trck.centers(3,:)',1,'r');
    plotmounts(Trck);
    
    % Include points from mounts-file in final point-cloud
    Terr.x = pts(:,1);
    Terr.y = pts(:,2);
    Terr.z = pts(:,3);
    
    Terr.g0 = g0; % approx. grid size (meters)
end

function D = readpointcloud(filename,pts,PRJ)
% [D,LOC] = READPOINTCLOUD(FILE,PTS,LOC,OPT)
% Read a FILE with three data columns X,Y,Z representing a set of scattered terrain points, with X,
% and Y in decimal degrees East/North (WGS84), and Z in mASL.
% Try to detect/fix small horizontal offsets between the points read in 
% the file, and a reference set of points PTS (expected to come from APPROXGROUNDMESH). 

    narginchk(3,3);
 
    D = readtxtfile(filename);
    D = cat(2,D.data{:});
    D = unique(D,'rows');
    D(any(~isfinite(D),2),:) = [];
    
    % Turn WGS84 to internal coordinate system
    [D(:,1),D(:,2)] = abs2prj(D(:,2),D(:,1),PRJ.origin(2),PRJ.origin(1),PRJ.rotation);
    D(:,3) = D(:,3) - PRJ.origin(3);
    
    k = convhull(D(:,1),D(:,2));
    assert(all(inpolygon(pts(:,1),pts(:,2),D(k,1),D(k,2))),'Points not in xyz region!');

    warning_disabler = naptime('terrainoffset:error','error'); %#ok<NASGU>
    
    % Check for errors in Mounts origin, ensure compatibility with XYZ cloud
    offset= findoffset(pts,D,{'mount positions','point-cloud data'});
    %assert(err(2) < warnerr,'Could not resolve differences between Mount-Positions and XYZ cloud.');

    D = D - offset;
end

function [pts,g0] = readgeotiff(filename,refpts,PRJ,opt)
% Read DEM from geotiff file, and check limits against mount-file coordinates
% Uses opt.rmin (with tolerance*) and opt.buffer

    % Min. RMS/MBE error to issue warning f(grid-lattice,refpts-median-dist)
    THRESH = getSimOption('terrain.tilethreshold');

    x0 = PRJ.origin(1);
    y0 = PRJ.origin(2);
    
    % Read file
    fprintf('\tOpening file...\n');
    Terr = geotiffread(filename);
    
    % Check boundaries
    if y0 > max(Terr.y) || y0 < min(Terr.y) || x0 > max(Terr.x) || x0 < min(Terr.x)
        error(['The provided geotiff file does not seem to cover the area of the project,',...
               'check the GEOTIFF file, and/or make sure that the Project''s origin in the',...
               '*.mounts file uses the syntax: #coords_center : LAT_N; LONG_E (decimal degrees).']);
    end

    % ... go even further, and check margins
    [x,y] = meshgrid(Terr.x([1,end]),Terr.y([1,end])); % get only boundaries, for now
    [x,y] = abs2prj(y,x,y0,x0,0);
    x = mean(x,1); y = flipud(mean(y,2)); % approx. edge midpoints

    margins = [x0 - x(1), x(2) - x0, y0 - y(1), y(2) - y0];
    msg = num2cell(margins/1000);
    msg = sprintf(['Geo-TIFF DEM extents [km]: %0.1f W, %0.1f E, %0.1f S, %0.1f N ',...
                   'from project circumcenter.'],msg{:});
    msg = {msg,sprintf('Recommended DEM radius: %0.1f km',opt.rmin/1000)}; 
    if any(margins < opt.rmin.*(1-THRESH))
        warning('%s\n%s',msg{:}); 
    else
        fprintf('\t%s\n\t%s\n',msg{:}); 
    end
    
    fprintf('\tTransforming to project-coordinates...\n');
    
    % Get all x,y pixel coordinates, and turn to internal coordinate system
    [x,y] = meshgrid(Terr.x,Terr.y);
    [x,y] = abs2prj(y,x,y0,x0,PRJ.rotation,'-approx');
    z = double(Terr.z);
    
    % Get lattice vectors
    u = cat(3,diff(x,1,2),diff(y,1,2));
    v = cat(3,diff(x,1,1),diff(y,1,1));
    g0(:,1) = squeeze(mean(mean(u,1),2));
    g0(:,2) = squeeze(mean(mean(v,1),2));
    
    % Resolve differences between the DEM and the elevation data in REFPTS
    % Use a local (CONVHULL+BUFFER) correction to achieve a smooth surface...

    fprintf('\tResolving differences with point-data...\n');
     
    offset = findoffset(refpts,[x(:),y(:),z(:)],{'ref. points','geotiff'});
    x = x - offset(1); y = y - offset(2); z = z - offset(3);

    % Find the region where we trust interpolations from scattered point-data
    % triangulation facets longer than s0 should not be trusted.
    s0 = mediandist(refpts(:,1),refpts(:,2));
    s0 = max(3*s0,norm(g0));
    p0 = polygon.outline(refpts(:,1),refpts(:,2),s0);  % ensures alpha-radius > s0, offsets s0
    p0 = polyout(p0,-s0/2);                            % bring down offset to s0/2
    
    % Create a larger buffer-polygon, that will be used to smoothly resolve differences
    p = polyout(p0,opt.buffer + 2*norm(g0),'r');
    near = inpolygon(x(:),y(:),p.x,p.y);
    inside = inpolygon(x(:),y(:),p0.x,p0.y);
    
    % Search the interpolated z values for refpts in the DEM
    IntZ = scatteredInterpolant(x(near),y(near),z(near),'natural');
    zt = IntZ(refpts(:,1),refpts(:,2));
        
    % Correct absolute offset of DEM data
    z = z - mean(zt) + mean(refpts(:,3));
    zt = zt - mean(zt) + mean(refpts(:,3));

    % Create an interpolant to smoothly resolve the differences between the DEM and mounts file
    IntZ = scatteredInterpolant([refpts(:,1);p.x'],[refpts(:,2);p.y'],[refpts(:,3) - zt;zeros(size(p.x))'],'natural');

    % Apply the interpolant-fix to any points inside the buffer polygon
    z(near) = z(near) + IntZ(x(near),y(near));
    
%     A0 = refinesquaregrid(ceil(norm(g0)/DMAX));
%     A = zeros([size(A0),numel(verynear)]);
%     for j = verynear'
%        A(:,:,j) = [x(j),y(j)] + A0*[u(:,j),v(:,j)]';
%     end

    % Include points from mounts-file in final point-cloud
    pts = cat(1,refpts,[x(~inside),y(~inside),z(~inside)]);
    g0 = norm(g0);
end

function [offset,err] = findoffset(pts,DEM,labels,varargin)
% [OFFSET,ERR] = FINDOFFSET(PTS,DEM,LABELS,...)
% Try to find an x,y offset that minimizes the RMS error between a set of points PTS and their
% corresponding interpolated values from a (larger) set of points DEM.
% The algorithm first evaluates a regular grid of offsets, then uses FMINSEARCH to refine the
% solution up to a given tolerance.
%
%   PTS - [n 3] array of points in x,y,z
%   DEM - [m 3] array of points, expected to be regularly spaced and cover an area larger than PTS
%   LABELS - 2-cellstring, describing the sets PTS and DEM
%
%   OFFSET - 3-vector [x,y,z] with the optimized result (must be added to PTS, to match DEM)
%   ERR - 2-vector [i,f] with initial and final RMS error.
%
% FINDOFFSET(..,'offsetradius',F = @(r)) search for best offset within F(d0).
% FINDOFFSET(..,'offsetlattice',N) use an N·N grid for the first (rough) search.
% FINDOFFSET(..,'tol',T) relative tolerance for fine-tuning (multiplied by d0)
% FINDOFFSET(..,'offsetP',P) consider offset as significant if the P-value of an F-test of the 
%     variances before and after the proposed correction is less than 1-P
% FINDOFFSET(..,'offsetwarning',T = @(g0)) issue a warning if the final STD error > T(g0)
%           
% (*) The DEM-grid lattice G0 is guessed by MEDIANDIST of a small subset of DEM, namely points
%     within d0 + 10·G0, where d0 is the plant 'diameter'.

    narginchk(3,9);

    opt = getSimOption('terrain'); % load complete terrain block (not all are used)
    opt.tol = getSimOption('RelTol');
    opt = getpairedoptions(varargin,opt,'restchk');
    
    % Get circumradius of plant, recenter everything around circumcenter 
    [d0,orig] = MinBoundCircle([pts(:,1),pts(:,2)]);    
    orig(3) = mean(pts(:,3));
    pts = pts - orig;
    DEM = DEM - orig;
    
    R0 = opt.offsetradius(d0); % search radius
    l0 = R0/opt.offsetlattice; % search lattice
    
    % Estimate the lattice distance of DEM grid
    r = hypot(DEM(:,1),DEM(:,2)); 
    near =  r < R0;
    g0 = mediandist(DEM(near,1),DEM(near,2));
    
    if ~opt.offsetcheck || all(abs(pts(:,3) - mean(pts(:,3))) < 0.001)
    % Offset-Check is actually disabled, just provide feedback
        
        offset = [0,0]; 
        near = r < d0 + g0;
        IntZ = scatteredInterpolant(DEM(near,1),DEM(near,2),DEM(near,3),'natural');
        err = std(pts(:,3)-IntZ(pts(:,1),pts(:,2)));
        
        % msg holds info for possible end warning
        if ~opt.offsetcheck, msg = 'Option terrain.offsetcheck is disabled. ';
        else,  msg = 'Points are flat to the millimeter, skipping terrain offset check. ';
        end
    else
    
        % Create a gridded-interpolant (for speed), note grid size of min(l0,g0)
        xq = linspace(-(R0+d0),R0+d0,ceil(2*(R0+d0)/min(l0,g0))); 
        [xq,yq] = ndgrid(xq,xq);
        near = r < R0 + d0 + g0;
        zq = griddata(DEM(near,1),DEM(near,2),DEM(near,3),xq,yq,'natural');
        IntZ = griddedInterpolant(xq,yq,zq,'linear');

        % Create a grid of point for first (rough) search
        xs = linspace(-(R0+d0),R0+d0,ceil(2*(R0+d0)/l0)); 
        [xs,ys] = ndgrid(xs,xs);
        incircle = xs.^2 + ys.^2 < R0^2;

        % Define the error function: f = std(z0 - Intz(x0+dx,y0+dy))im
        % meandev = @(x) mean(abs(x-mean(x)));
        % stderrf = @(x) meandev(pts(:,3)-IntZ(pts(:,1)+x(1),pts(:,2)+x(2)));
        stderrf = @(x) std(pts(:,3)-IntZ(pts(:,1)+x(1),pts(:,2)+x(2)));
        err = stderrf([0,0]);

        if R0 == 0 || err == 0
        % No need to search
            offset = [0,0];
            msg = 'Zero search radius!, skipping terrain offset check. ';
        else   
            % Evaluate errors when applying offset at grid-points xs, ys in circle
            c = NaN(size(xs));
            c(incircle) = arrayfun(@(x,y) stderrf([x,y]),xs(incircle),ys(incircle));

            % Find local minimum near (0,0) offset, and add it to the list of possible values
            [offset,c0] = fminsearch(stderrf,[0,0],struct('TolX',d0*opt.tol));
            
            % Call FMINSEARCH using the best value as seed (it could well be 0,0)
            [~,ic] = min([c(:);c0],[],'omitnan');
            if ic < numel(c)
                offset = fminsearch(stderrf,[xs(ic),ys(ic)],struct('TolX',d0*opt.tol));
            end

            e0 = IntZ(pts(:,1),pts(:,2))-pts(:,3);
            e2 = IntZ(pts(:,1)+offset(1),pts(:,2)+offset(2))-pts(:,3);
            err(2) = std(e2); % min. found error

            % Get the probability that errors e0 and e2 have different variances
            P = 1 - vartestn([e0,e2],'TestType','LeveneAbsolute','Display','off');
            
            % fraction of points for which err(xq,yq) < err(0,0) 
            locality = 1-nnz(c < err(1))/nnz(incircle);

            if norm(offset) > d0*opt.tol && P > opt.offsetP && runningfromUI()
            % Ask user whether to apply significant, non-zero offset
                if runningfromUI()
                    GUIfigure('terrain_offset','Terrain Offset Check','2:1'); clf(); 
                    subplot(1,2,2); hold on;
                    c(~incircle) = NaN;
                    contour(xs,ys,c,linspace(err(2),min(2*err(1)-err(2),max(c(:),[],'omitnan')),10));
                    colorbar();
                    axis equal; grid on;
                    axis([-1 1 -1 1]*ceil(R0/10)*10);
                    plot(0,0,'k+','markersize',20);
                    plot(offset(1),offset(2),'r+','markersize',20);
                    title(sprintf('STD (%s - %s)',labels{1},labels{2}))
                    xlabel('X offset (m)');
                    ylabel('Y offset (m)');
                    
                    subplot(1,2,1); hold on;
                    axis equal; grid on;
                    axis([-1 1 -1 1]*ceil(R0/10)*10);
                    DT = delaunayTriangulation(DEM(near,1:2));
                    trisurf(DT.ConnectivityList,DEM(near,1),DEM(near,2),DEM(near,3)-mean(e0),'edgealpha',0.05);
                    plot3(pts(:,1),pts(:,2),pts(:,3),'r.');
                    title('Current position');
                    xlabel('X (m)');
                    ylabel('Y (m)');
                end
                msg = sprintf(['An offset of %0.1f m X, %0.1f m Y would reduce the error ',...
                  'between %s and %s from %0.2f m to %0.2f m (%0.1f%% significance). ',...
                  'The current position is better than %0.2f%% of locations within %0.1fkm. ',...
                  'Apply the offset?'], offset(1),offset(2),labels{1},labels{2},...
                   err(1),err(2),P*100,locality*100,R0/1000);
                switch optquestdlg(msg,'Terrain Offset','Apply','Ignore','Ignore')
                    case 'Apply'
                        msg = sprintf('Offset of %0.1f m X, %0.1f m applied to %s. ',...
                            offset(1),offset(2),labels{2});
                        err = err(2);
                    case 'Ignore'
                        msg = sprintf('Detected %0.1f m X, %0.1f m offset ignored. ',...
                            offset(1),offset(2));
                        offset = [0,0]; 
                        err = err(1);
                    otherwise
                        error('Interrupted by user');
                end
            else
                msg = sprintf(['Failed to find non-zero offset within search radius (%0.1fkm) ',...
                    'that would result in a significant fit improvement. ',...
                    'The current position is better than %0.1f%% of neighboring locations. ',...
                    'Best guess was %0.1f m X, %0.1f m Y (STD = %0.2f m, P = %0.1f%%). '],...
                    R0/1000,locality*100,offset(1),offset(2),err(2),P*100);
                err = err(1);
                offset = [0,0]; 
            end
        end
    end
    offset(3) = mean(IntZ(pts(:,1)+offset(1),pts(:,2)+offset(2))-pts(:,3));
    
    % Issue warning only if error is above given threshold
    if err > opt.offsetwarning(g0)
        warning('terrainoffset:error', [msg,...
             'Differences between %s and %s (STD = %0.2f m) persist. Try resolving ',...
             'them using GIS/3D-modeling software.'],labels{1},labels{2},err);
    end
end

function [h,v] = mediandist(x,y,z)
% h = MEDIANDIST(x,y) gives a representative distance metric for the set of scattered points x,y, 
% as the median of the edge-lengths of their (flat) Delaunay triangulation.
% [h,v] = MEDIANDIST(x,y,z) - additionally provides the median absolute z-distance for all edges.

    DT = delaunayTriangulation(x,y);
    E = edges(DT);
    d = rssq(DT.Points(E(:,1)) - DT.Points(E(:,2)),2);
    h = median(d);
    
    if nargin > 2 && nargout > 1
        assert(numel(z) == size(DT.Points,1),'Mean-vertical distance requires unique points');
        v = median(abs(z(E(:,1)) - z(E(:,2))));
    end
end

function varargout = plotmounts(Trck)
% Minimal Layout-plot (only mount patches, sun at zenith)

    [F,V] = plantlayout(Trck,0,90-Trck.origin(2),'eq0');
    h = patch('Faces',F,'Vertices',V,'FaceColor',[0.3,0.3,0.8],'edgecolor','w','edgealpha',0.1);
    
    if nargout > 0, varargout{1} = h; end
end