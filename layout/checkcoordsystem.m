function [x,y,z,S,info] = checkcoordsystem(x,y,z,varargin)
% [X,Y,Z,S] = CHECKCOORDSYSTEM(LON,LAT,ALT,[P],..,['output','prj'])
% [LAT,LON,ALT,S] = CHECKCOORDSYSTEM(X,Y,Z,[P],..,'output','abs')
%   Parse coordinate-system specification P (alternatively provided/modified through name,value
%   pairs), and transform input coordinates from system P onto project-centered (rotated) Transverse
%   Mercator ('prj') coordinates or WGS84 decimal degree ('abs' coordinates).
%   Return also a structure S with the final projection system specs.
%
%   Valid fields for P / S are:
%
%   system: 
%       'custom'/'prj': (default) local projection (assumed conformal & equidistant), centered on
%           a given ORIGIN and with a given ROTATION from geographical North (see below).
%       'WGS84' / 'EPSG:4326': x,y coordinates read as decimal degrees longitude & latitude. 
%           Make sure to provide enough decimals (7-8 for cm/mm precision).
%       'UTM XXZ': x,y coordinates in Universal Transverse Mercator projection (zone XXZ, e.g. 30T)
%       '+proj=..': custom PRJ.4 projection string (requires installed OSGeo/proj). See PROJ.
%
%   latitude, longitude, altitude: alternatives to provide origin (below). 
%
%   origin: 3 vector [x0,y0,z0] / [lon0,lat0,z0] or string 'lat0 N, lon0 E, z0]'. 
%       For CUSTOM coords_system, lat0,l on0 must be WGS84 decimal-degree latitude (N) and 
%       longitude (E) (at least 6-decimal precision) of the project-coordinate-system-origin. 
%       For UTM projections, x0,y0 can be a common offset (meters) to be added to x,y coordinates
%       before conversion to WGS84.
%       Optional z0 will be read as altitude in mASL, and simply added to z.
%
%   rotation: decimal-degree rotation of project-coordinate-system CCW from geographical North.
%       rotation = 0 <-> Y-axis = North, X-axis = East; rotation = 90 <-> Y-axis = West
%
% ..,'-recenter' - check that X,Y,Z coordinates are reasonably* centered around the specified 
%       origin, and adjust the origin otherwise.
%
% ..,'-plot','-inland' - passed to PARSELOCATION
%
% See also: IMPORTPLANTLAYOUT, PROJ, ABS2PRJ, PRJ2ABS, PARSELOCATION

    ALIAS = {'rotation',{'rotation','coords_rotation'};
             'origin',{'center','origin','coords_center'};
             'system',{'system','coords_system','proj','projection'};
             'lat',{'lat','latitude'};
             'lon',{'lon','longitude'};
             'alt',{'alt','altitude'}};
         
    nout = nargout; info = {}; % used to capture warnings if nargout > 4
         
    [opt,varargin] = getflagoptions(varargin,{'-inland','-recenter','-plot'});     
    opt.output = 'prj';
    opt.offset_tolerance = 1.0;    % Allowed distance (in units of mean mount distance from minimax
                                   % center, before recentering.
    [opt,varargin] = getpairedoptions(varargin,opt);
    opt.output = parselist(opt.output,{'prj',{'prj','project','local'};
                                       'abs',{'abs','absolute','wgs84','deg'}},'output type');
    
    [params,varargin] = getpairedoptions(varargin,cat(2,ALIAS{:,2}));
    if ~isempty(varargin)
       assert(isscalar(varargin) && isstruct(varargin{1}),'Unrecognized arguments');
       params = completestruct(params,varargin{1});
    end
    params = renamefields(params,fliplr(ALIAS),'-ignorecase');
    
    DEF = struct('system','','rotation',[],'origin',[]);
    if all(isfield(params,{'lat','lon'})), DEF.origin = [params.lon,params.lat]; end
    params = completestruct(params,DEF);

    system = params.system;
    if ~isempty(system)
        assert(isempty(system) || ischar(system),'Expecting char system');
    end
    
    rotation = params.rotation;
    if ~isempty(rotation)
        if ischar(rotation), rotation = str2double(rotation); end
        validateattributes(rotation,{'numeric'},{'real','scalar','finite'},'','rotation');
    end
    
    origin = params.origin;
    if ~isempty(origin)
        if ischar(origin), origin = str2origin(origin);
        else
            validateattributes(origin,{'numeric'},{'real','vector','finite','nonempty'},'','origin')
            assert(any(numel(origin) == [2,3]),'Expecting 2-3 element origin vector');
        end
        if numel(origin) < 3
            if isfield(params,'alt') , origin(3) = params.alt;
            else, origin(3) = NaN;
            end
        end
    end
    
    if isempty(system) || any(strcmpi(system,{'prj','project','local','custom','TM'}))
    % Input coords are custom TM projection centered on origin
    
        assert(~isempty(origin) && abs(origin(2)) <= 90 && abs(origin(1)) <= 180,...
            'Invalid ORIGIN. Make sure to use LAT_N; LONG_E (decimal degrees)');
    
        lon0 = origin(1);
        lat0 = origin(2); 
        
        S.origin = origin;
        if isempty(rotation), rotation = (lat0 > 0)*180; end

        switch opt.output
        case 'prj'
            S.system = 'custom';
            S.origin = origin;
            S.rotation = rotation;
            % x = x, y = y, z = z
        case 'abs'
            S.system = 'WGS84';
            [y,x] = prj2abs(x,y,lat0,lon0,rotation);
            if ~isnan(origin(3)), z = z + origin(3); end
        end
    else
        if regexpi(system,'^(epsg.{0,3}4326|wgs84|abs)$')  
        % WGS84 (decimal degree) projection
        
            assert(all(abs(y) <= 90) && all(abs(x) <= 180),...
                'Expecting decimal degree coordinates x (longitude) and y (latitude)');

            lon = x;
            lat = y;
            % z = z;

        elseif regexpi(system,'utm.?\d{1,2}[a-z]{1}') 
        % UTM projection
        
            zone = regexprep(upper(system),'UTM[^\d]?(\d{1,2}[A-Z]{1})','$1');
            [lat,lon] = utm2deg(x+origin(1),y+origin(2),zone);
            
        else
            % try
                [lon,lat] = proj(system,x,y);              
            % catch ERR
            %     error('Coordinate System (%s) is not recognized!',system)
            % end
        end 
        
        lon0 = round(mean(lon,'omitnan'),6);
        lat0 = round(mean(lat,'omitnan'),6);
    
        switch opt.output
        case 'prj'
        % Custom TM projection centered on (rounded) mean point
        
            if isempty(rotation), rotation = (lat0 > 0)*180; end
         
            S.system = 'custom';
            S.rotation = rotation;
            S.origin = [lon0,lat0,origin(3)];
            [x,y] = abs2prj(lat,lon,lat0,lon0,rotation);
            if ~isnan(origin(3)), z = z - origin(3); end
        case 'abs'
            if ~isempty(rotation) && rotation ~= 0
                warningOrInfo('norotation''ignoring rotation = %0.0f°',rotation);
            end
            S.system = 'WGS84';
        end
    end
            
    if opt.inland || opt.plot
        if opt.inland, args = {'-inland'}; else, args = {}; end
        if opt.plot, args{end+1} = '-plot'; end
        parselocation(struct('latitude',lat0,'longitude',lon0,'altitude',0),args{:});
    end
     
    if ~opt.recenter || ~strcmp(opt.output,'prj'), return; end
    
    % Check if the coordinate system's origin is reasonably centered...    
    c = [(min(x)+max(x))/2,(min(y)+max(y))/2];
    d0 = hypot(max(x)-min(x),max(y)-min(y));

    if norm(c) > d0*opt.offset_tolerance
    % ... if it's not, offset by circumcenter (round-up to 1e-6 decimal degrees)
    
        [lat0,lon0] = prj2abs(c(1),c(2),lat0,lon0,rotation);
        lon0 = round(lon0,6);
        lat0 = round(lat0,6);
                
        warningOrInfo('recenter',['Re-centering coordinate system, ',...
            'from: %0.6f°N, %0.6f°E to %0.6f°N, %0.6f°E'],S.origin(2), S.origin(1),lat0,lon0);
   
        originoffset = abs2prj(S.origin(2),S.origin(1),lat0,lon0,rotation);
        S.origin(1:2) = [lon0,lat0];
        x = x - originoffset(1);
        y = y - originoffset(2);
    end
    
    if isnan(origin(3)), return; end
    c = min(z);
    d0 = max(z)-min(z);

    if abs(c) > d0*opt.offset_tolerance
        c = round(c,1);
        warningOrInfo('height','Applying altitude offset: from %0.1f to %0.1f mASL',S.origin(3),c);
        z = z - c; 
        S.origin(3) = S.origin(3) + c;
    end
    
    function warningOrInfo(id,msg,varargin)
        if nout > 4
            info{end+1} = ['warning: ' sprintf(msg,varargin{:})];
        else
            warning(['checkcoordsystem:' id],msg,varargin{:});
        end
    end
end

function c = str2origin(s)
    s = upper(s);
    try
        if contains(s,{'N','S','E','W'})
            s = strsplit(s,';');
            if isscalar(s), s = strsplit(s,','); end
            if isscalar(s), s = strsplit(s,' '); end
            assert(any(numel(s) == [2,3]),'Expecting 2-3 vector');

            if contains(s{1},{'N','S'}) && contains(s{2},{'E','W'}), s(1:2) = s([2,1]); end
            c = cellfun(@dms2deg,s(1:2));
            if numel(s) > 2
                s{3} = regexprep(s{3},'^\s*(\d+)\s*.*','$1');
                c(3) = str2double(s{3}); 
            end
        else
            c = str2num(s); %#ok<ST2NM>
            assert(any(numel(c) == [2,3]),'Expecting 2-3 vector');
            c(1:2) = c([2,1]);
        end
        c = c(:)';
        validateattributes(c(1:2),{'numeric'},{'real','finite'})
    catch ERR
       error('Failed to understand origin string: %s',ERR.message); 
    end
end