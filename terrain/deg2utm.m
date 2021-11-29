function  [x,y,zone,letter] = deg2utm(Lat,Lon,varargin)
% [X,Y,Z,L] = DEG2UTM(LAT,LON) convert lat/lon arrays into UTM coordinates (WGS84).
% [X,Y,Z,L] = DEG2UTM(LAT,LON,ZONE) force the use of a given ZONE.
% [X,Y,Z,L] = DEG2UTM(...,'name',Value) modify custom parameters individually.
%
% [X,Y,Z,L] = DEG2UTM(LAT,LON,'-multizone') - allow points to use different (automatic) UTM zones,
%   otherwise the MODE of the automatic zone assignment will be used.
% [LAT,LON] = DEG2UTM(...,'-approx') - use approximate projection (~ 2 cm accuracy) [2]
%
%   LAT: Latitude array.   Degrees North
%   LON: Longitude array.  Degrees East. NOTE: size(LAT) must equal size(LON).
%   ZONE: 1x3 or nx3 array of zone codes (see '-multizone' option below) with format ##A, eg. '30T'.
%
%   X,Y: size(LAT) arrays of UTM coordinates.
%   Z: UTM zone number(s)
%   L: UTM zone letter(s)
%
%   'name',Value parameters (custom Transverse Mercator projection):
%
%       lat_0 - Latitude of origin
%       lon_0 - Central zone meridian (required if zone is not provided)
%       x_0 - False easting, default is 500000
%       y_0 - False northing, default is 1e7·(lat_0 < 0)
%       ellps - Ellipsoid, default is WGS84
%       k_0 - Scaling factor, default 0.9996
%       units - Default 'm'
% 
% [1] PROJ contributors (2020). PROJ coordinate transformation software library. 
%     Open Source Geospatial Foundation. URL https://proj.org/.
% [2] Modified from: Rafael Palacios, Aug/06. Universidad Pontificia Comillas, Madrid.
% 
% See also: UTM2DEG, PROJ, PRJ2ABS

    narginchk(2,16);
    validateattributes(Lat,{'numeric'},{'real','>=',-90,'<=',90},'','Lat');
    validateattributes(Lon,{'numeric'},{'real','>=',-180,'<=',180},'','Lon');
    [Lat,Lon] = compatiblesize(Lat,Lon); 
    
    [opt,varargin] = getflagoptions(varargin,{'-multizone','-approx'});
    opt.ellps = 'WGS84';
    opt.lat_0 = 0;
    opt.lon_0 = [];
    opt.x_0 = 500000;
    opt.y_0 = 0;
    opt.k_0 = 0.9996;
    opt.units = 'm';
    [opt,varargin,isdef] = getpairedoptions(varargin,opt);
    cellfun(@(f) validateattributes(opt.(f),{'char','string'},{'nonempty'},'',f),{'ellps','units'});

    if isempty(varargin)
        zone = ''; 
    else
        assert(isscalar(varargin),'Unrecognized arguments');
        zone = varargin{1}; 
    end
    
    % Check for OSGeo/proj on system
    systemproj = ~opt.approx;
    if systemproj
        try proj('+proj=ortho',0,0); catch, systemproj = false; end
    end
    
    if nargout > 3
        letter = zoneletter(Lat,Lon);
        if ~opt.multizone, letter = mode(letter); end
    end
    
    trueUTM = all(cellfun(@(f) isdef.(f),setdiff(fieldnames(opt),'ellps')));
        
    south = Lat < 0;
    if ~isempty(zone)
    % Provided zone (i.e. char-code zone)
        
        if ischar(zone)
            assert(size(zone,2) == 2 || ...
                ( size(zone,2) == 3 && all(ismember(upper(zone(:,3)),'A':'Z')) ),...
                'Wrong UTM-ZONE format, expecting form ##[A]');

            if size(unique(zone,'rows'),1) > 1, opt.multizone = true; end

            if size(zone,2) > 2
                letter = zone(:,3);
                south = upper(letter) < 'N';
            end
            zone = str2num(zone(:,1:2)); %#ok<ST2NM>
        else
            assert(isnumeric(zone) && isreal(zone),'Expecting char or numeric UTM zone');
        end
        assert(~any(zone < 1 | zone > 60 | mod(zone,1) ~= 0),'Invalid UTM-ZONE number');
        
        if ~trueUTM
            warning('UTM zone(s) overridden by custom TM parameters, not a true UTM!')
        end
        
        z = zonenumber(Lon,Lat);
        if any(rem(abs(zone - z),59) > 1 & Lat > -80 & Lat < 84)
        % allow overlap to neighboring zones
            warning('Enforced UTM-ZONE(s) seem too far from some points, double-check your pick');
        end
        
        if isdef.lon_0, opt.lon_0 = ((zone * 6) - 183); end
        if isdef.y_0, opt.y_0 = south*10000000; isdef.y_0 = false; end
    
    elseif isdef.lon_0
    % Custom projection parameters

        [zone,opt.lon_0] = zonenumber(Lon,Lat);
    end
            
    if trueUTM && systemproj
    % proj +proj=utm +zone=XX ...

        if ~opt.multizone
            zone = mode(zone(:)); 
            south(:) = mode(south(:));
        else
            [south,zone] = compatiblesize(south,zone,Lat);
        end
        
        prj4 = arrayfun(@(z) sprintf('+proj=utm +ellps=%s +zone=%d',opt.ellps,z),zone,'unif',0);
        prj4(south) = cellfun(@(s) [s ' +south'],prj4(south),'unif',0);

        [x,y] = proj(prj4,Lon,Lat);
        return;
    end

    if isdef.y_0
        south = Lat < 0;
        opt.y_0 = 1e7*south; 
    end
    
    parsestruct(opt,{'lat_0','lon_0','x_0','y_0','k_0'},'-n','-r','nonempty');
    validateattributes(opt.lat_0,{'numeric'},{'>=',-90,'<=',90},'','lat_0');

    [lat_0,lon_0,x_0,y_0,k_0] = deal(opt.lat_0,opt.lon_0,opt.x_0,opt.y_0,opt.k_0);
    [lat_0,lon_0,x_0,y_0,k_0] = compatiblesize(lat_0,lon_0,x_0,y_0,k_0);
    
    if ~opt.multizone && ~isscalar(lat_0)
        x = num2cell(mode([lat_0,lon_0,x_0,y_0,k_0],1)); 
        [lat_0,lon_0,x_0,y_0,k_0] = deal(x{:});
    end

    if systemproj
    % Use OSGeo/proj, when available. 
    
        prj4 = [ '+proj=tmerc +ellps=' opt.ellps ' +units' opt.units ];
        prj4 = arrayfun(@(varargin) sprintf(...
            [prj4 ' +lat_0=%0.8f +lon_0=%0.8f +x_0=%0.3f +y_0=%0.3f +k_0=%0.6f'],varargin{:} ),...
            lat_0,lon_0,x_0,y_0,k_0,'unif',0);
        
        [x,y] = proj(prj4,Lon,Lat);
        return;     
    end

    if ~strcmpi(opt.ellps,'WGS84')
       warning('Custom ellipsoids require PROJ system package: using WGS84');
    end
    
    % The following code was modified from: 
    % Rafael Palacios (2021). utm2deg, MATLAB Central File Exchange.
    % <https://www.mathworks.com/matlabcentral/fileexchange/10914-utm2deg>

    sa = 6378137.000000 ; sb = 6356752.314245;
    % e2 = (sa^2 - sb^2)/sb^2;
    % c = (sa^2) / sb;
    e2 = (sa/sb)^2 - 1;
    e = sqrt(e2);
    c = (sa/sb)*sa;

    LatR = Lat * (pi / 180);
    LonR = Lon * (pi / 180);
    
    deltaS = LonR - (lon_0 * (pi / 180));
    
    maxdS = max(abs(deltaS))*180/pi;
    if maxdS > 6
        warning([ 'Distances from central meridian of up to %0.1f°,',...
            ' verify your choice of zone, or consider using -multizone'],maxdS); 
    end
    
    [y,x] = UTM_prj(deltaS,LatR);
    
    % Add false easting/northing
    if lat_0 ~= 0
        y_0 = y_0 - UTM_prj(0,lat_0*pi/180);
    end
    %y(insouth) = 9999999+y(insouth);
    x = x + x_0;
    y = y + y_0;
    
    function [y,x] = UTM_prj(deltaS,LatR)
        a = cos(LatR) .* sin(deltaS);
        epsilon = 0.5 * log((1 + a) ./ (1 - a));
        nu = atan(tan(LatR) ./ cos(deltaS)) - LatR;
        % v =  c ./ sqrt(1 + e2 * cos(LatR).^ 2) * 0.9996;
        v = c./hypot(1,e*cos(LatR))*0.9996;
        ta = (e2 / 2) * epsilon.^2 .* cos(LatR).^2;
        a1 = sin(2 * LatR);
        a2 = a1 .* cos(LatR).^2;
        j2 = LatR + (a1 / 2);
        j4 = (3 * j2 + a2) / 4;
        j6 = (5 * j4 + a2 .* cos(LatR).^2) / 3;
        alfa = (3 / 4) * e2;
        beta = (5 / 3) * alfa.^2;
        gama = (35 / 27) * alfa.^3;
        Bm = 0.9996 * c .* (LatR - alfa .* j2 + beta .* j4 - gama .* j6);
        x = epsilon .* v .* (1 + ta / 3);
        y = nu .* v .* (1 + ta) + Bm;
    end
end

function [n,lon_0] = zonenumber(Lon,Lat)

    n = discretize(Lon,-180:6:180);
    %n(~isfinite(n)) = 1;
    
    % exceptions around Norway
    exc = Lat >= 56 & Lat < 64;
    if any(exc)  
        n(exc & Lon >= 3 & Lon < 12) = 32;
    end
    exc = Lat >= 72 & Lat < 84;
    if any(exc)
        exc = exc & Lon >= 0 & Lon < 42;
        n(exc) = discretize(Lon(exc),[0,9,21,33,42])*2+29;
    end
    
    lon_0 = (n * 6) - 183;
end

function L = zoneletter(Lat,Lon)
% Get standard latitudinal zone letter

    ZoneLetterList = ('ABCDEFGHJKLMNPQRSTUVWXYZ#')';
    
    L = discretize(Lat,[-90,-80:8:72,84,90])+1;
    
    % Regions A, B, Y, Z near poles
    L(L == 2 & Lon < 0) = 1;
    L(L == 23 & Lon > 0) = 24;
    
    L(~isfinite(L)) = 25; % Flag NaNs as #
    
    L = ZoneLetterList(L);
end


