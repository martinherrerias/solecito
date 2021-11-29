function varargout = getgeotiff(lat,lon,varargin)
% F = GETGEOTIFF(LAT,LON,RADIUS)
% Attempt to find (retrieve from URL if necessary)* the required ASTER 1×1° DEM tiles that  
% cover a circle of a given RADIUS from LAT, LON; then merge them and clip them into a single
% geotiff file.
%
% F = DEMRASTERTILES(..,'thresh',F) tiles that cover less than a fraction F of the circle will
%   not be included in the list (default is 5%).
%
% F = DEMRASTERTILES(..,'local') search for required tiles locally, do not open any URL windows.
%
% (*)FUTURE: find an API for auto-download. The best we can do so far is open a link:
% https://search.earthdata.nasa.gov/search/?q=ASTER&circle[0]=...
    
    [opt,remargs] = getflagoptions(varargin,{'local'});
    opt.thresh = getSimOption('terrain.tilethreshold'); 
    [opt,remargs] = getpairedoptions(remargs,opt);
    if ~isempty(remargs), R = remargs{1}; else, R = getSimOption('terrain.demradius'); end
    
    [tilekeys,url,c] = demrastertiles(lat,lon,R,opt.thresh);

    findfiles = @(K) cellfun(@(k) dir(['*' k '*.tif']),K,'unif',0);

    % Check if files for those tiles are already available
    files = findfiles(tilekeys);
    missing = cellfun(@isempty,files);
    
    if any(missing) && runningfromUI()
        if ~opt.local
            % Open browser download links
            web(url,'-browser');
            % Keep nagging the user until every file is available
            msg = ['Please download the required ASTER-V2-DEM tiles. Make sure to ',...
                'copy all *.tif files in the current directory, and include the ',...
                shortliststr(tilekeys(missing),'tile-key'),...
                ' as part of the file-name(s). Close this dialog when you''re ready to continue'];
            switch questdlg(msg,'Missing DEM tiles','Continue','Cancel','Continue')
            case {'Cancel',''}
                   error('getgeotiff:cancelled','Could not retrieve all the required *.tif files'); 
            end
        end

        % Keep launching UI dialogs to select missing files until the user cancels
        msg = ['Select missing ' shortliststr(tilekeys(missing),'tile','colon',':')];
        ff = pickfile('*.tif',Inf,msg,'ui',1);
        if ~isempty(isempty(ff))
            for j = find(missing)
                useful = contains(ff,tilekeys(j));
                if any(useful)
                   files = cat(1,files,ff(useful));
                   missing(j) = false;
                end
            end
            files = uniquecell(files);
        end
    end
    if any(missing) && ~all(missing)
        msg = sprintf('%d of the %d recommended tiles are missing, merge anyway?',nnz(missing),numel(tilekeys));
        switch questdlg(msg,'Missing DEM tiles','Yes','Cancel','Cancel')
        case {'Cancel',''}
               error('getgeotiff:cancelled','Could not retrieve all the required *.tif files');
        end
    end
    files = files(~missing);
    
    g = cellfun(@(f) geotiffread(fullfile(f.folder,f.name)),files);
    
    % Merge all tiles into a single image G
    [x,~,ix] = unique([g.x]);
    [y,~,iy] = unique([g.y]);
    ix = reshape(ix,[],numel(g));
    iy = reshape(iy,[],numel(g));
    G = nan(numel(y),numel(x));
    for j = 1:numel(g), G(iy(:,j),ix(:,j)) = g(j).z; end
    
    % Clip to circle boundaries
    c = offsetpolygon(c,1/3600);
    inxbounds = x > min(c.x) & x < max(c.x);
    inybounds = y > min(c.y) & y < max(c.y);
    x = x(inxbounds);
    y = y(inybounds);
    G = G(inybounds,inxbounds);
    
    % Generate a file-name of the form: 'ASTGTM2_21N20_109W14_R50km_dem.tif'
    deg2dm = @(x,k) sprintf('%d%c%d',floor(abs(x)),k,round(60*rem(abs(x),1)));  
    filename = [deg2dm(lat,(lat < 0)*('S' - 'N')+'N') '_' deg2dm(lon,(lon < 0)*('W' - 'E')+'E')];
    filename = sprintf('ASTGTM2_%s_R%dkm_dem.tif',filename,round(R/1000));
    
    % Write clipped geotiff
    G = flipud(G);
    bbox = [min(x),min(y);max(x),max(y)];
    geotiffwrite(filename,bbox,G,16,struct('GTRasterTypeGeoKey',2));
    fprintf('\tMerged Geo-TIFF written to: %s\n',relativepath(filename));
    if nargout > 0, varargout{1} = filename; end
end

function [K,URL,c] = demrastertiles(lat,lon,radius,thresh)
% [K,U] = DEMRASTERTILES(LAT,LON,RADIUS,THRESH)
% Returns a cellstr K of tile-keys (e.g. 'N52E012') and a cellstr U of URL-addresses to download
% the required ASTER V2 1×1° tiles that cover a circle of a given RADIUS from LAT, LON.
% Tiles that cover less than a fraction THRESH of the circle will not be included in the list
% (default is 5%).
    
    URL = 'https://search.earthdata.nasa.gov/search?q=ASTER';
    URL = [URL, sprintf('&circle[0]=%0.6f%%2C%0.6f%%2C%d',lon,lat,radius)];
    URL = [URL, sprintf('&m=%f!%f!',lat,lon),'7!1!0!0%2C2'];
    
    % Custom (project-centered) Transverse-Mercator Projection
    % UTM_PRJ = round([lat,lon,0,0],3);
    % [x0,y0,UTM_PRJ] = deg2utm(lat,lon,UTM_PRJ);
    x0 = 0; y0 = 0;

    % Get 1° tiles required for RADIUS coverage:
    c = polygon(360);
    c.x = c.x*radius + x0; 
    c.y = c.y*radius + y0;
    % [c.y,c.x] = utm2deg(c.x,c.y,UTM_PRJ); 
    [c.y,c.x] = prj2abs(c.x,c.y,lat,lon,0);
    
    [x,y] = meshgrid(floor(min(c.x)):floor(max(c.x)),floor(min(c.y)):floor(max(c.y)));
    
    for j = numel(x):-1:1
       q = intersectpolygons(c,polygon([0,1]+x(j),[0,1]+y(j)));
       if q.area < c.area*thresh, continue; end
       K{j} = sprintf('%c%02d%c%03d',(y(j) < 0)*('S' - 'N')+'N',abs(y(j)),(x(j) < 0)*('W' - 'E')+'E',abs(x(j)));
    end   
    K(cellfun(@isempty,K)) = [];
    K = K(:);
end
