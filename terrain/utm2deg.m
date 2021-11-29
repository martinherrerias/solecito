function  [Lat,Lon] = utm2deg(x,y,varargin)
% [LAT,LON] = UTM2DEG(X,Y,ZONE) convert UTM arrays into lat/lon coordinates (WGS84).
% [LAT,LON] = UTM2DEG(X,Y,N,SOUTH) provide zone as number, and boolean SOUTH (Lat < 0)
% [LAT,LON] = UTM2DEG(...,'name',Value) modify custom parameters individually.
% [LAT,LON] = UTM2DEG(...,'-approx') - uses approximate projection (~ 2 cm accuracy) [2]
%
%   X,Y: equal-size arrays of UTM coordinates.
%   ZONE: 1x3 string (single zone) or nx3 char-array (one zone for each X,Y), where each row has
%       the form '##A', eg. '30T'.
%   LAT: Latitude array.   Degrees North
%   LON: Longitude array.  Degrees East
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
% See also: DEG2UTM, PROJ, ABS2PRJ

    narginchk(3,16);
        
    [opt,varargin] = getflagoptions(varargin,{'-approx'});
    opt.ellps = 'WGS84';
    opt.lat_0 = 0;
    opt.lon_0 = [];
    opt.x_0 = 500000;
    opt.y_0 = 0;
    opt.k_0 = 0.9996;
    opt.units = 'm';
    
    [opt,varargin,isdef] = getpairedoptions(varargin,opt);
    cellfun(@(f) validateattributes(opt.(f),{'char','string'},{'nonempty'},'',f),{'ellps','units'});
    
    % Check for OSGeo/proj on system
    systemproj = ~opt.approx;
    if systemproj
        try proj('+proj=ortho',0,0); catch, systemproj = false; end
    end
    
    trueUTM = all(cellfun(@(f) isdef.(f),setdiff(fieldnames(opt),'ellps')));
    south = [];
            
    if ~isempty(varargin)
    % Specified Zone (UTM)
    
        assert(numel(varargin) <= 2,'Unexpected arguments');
        if numel(varargin) > 1, south = varargin{2}; end
        
        if ~trueUTM
            warning('UTM zone(s) overridden by custom TM parameters, not a true UTM!')
        end
 
        if ischar(varargin{1})
        % Standard (i.e. char-code zone)
        
            zone = strtrim(varargin{1}); clear varargin

            assert(size(zone,2) == 2 || ...
                ( size(zone,2) == 3 && all(ismember(upper(zone(:,3)),'A':'Z')) ),...
                'Wrong UTM-ZONE format, expecting form ##[A]');
        
            if isempty(south) && size(zone,2) > 2, south = zone(:,3); end
            zone = str2num(zone(:,1:2)); %#ok<ST2NM>
        else
            zone = varargin{1};
            assert(isnumeric(zone) & isreal(zone),'Expecting char or numeric ZONE');
        end
            
        if ~isempty(south)
            if ischar(south), south = upper(south) < 'N'; end
            compatiblesize(south,x);
            assert(all(south == 1 | south == 0),'Expecting boolean array SOUTH');
        end
        
        % assert(~any(zone < 1 | zone > 60 | mod(zone,1) ~= 0),'Invalid UTM-ZONE number');

        if isdef.lon_0, opt.lon_0 = ((zone * 6) - 183); end
    end
    
    if isdef.y_0
        assert(~isempty(south),'Zone letter, custom y_0, or boolean SOUTH required');
        opt.y_0 = south*10000000; 
    end
    
    if trueUTM && systemproj
    % proj -I +proj=utm +zone=XX ...
        
        south = compatiblesize(south,x);
        
        prj4 = arrayfun(@(z) sprintf('-I +proj=utm +ellps=%s +zone=%d',opt.ellps,z),zone,'unif',0);
        prj4(south) = cellfun(@(s) [s ' +south'],prj4(south),'unif',0);
        
        [Lon,Lat] = proj(prj4,x,y);
        return;
    end
    
    parsestruct(opt,{'lat_0','lon_0','x_0','y_0','k_0'},'-n','-r',@(x) ~isempty(x));
    validateattributes(opt.lat_0,{'numeric'},{'>=',-90,'<=',90},'','lat_0');

    [lat_0,lon_0,x_0,y_0,k_0] = deal(opt.lat_0,opt.lon_0,opt.x_0,opt.y_0,opt.k_0);
    [lat_0,lon_0,x_0,y_0,k_0] = compatiblesize(lat_0,lon_0,x_0,y_0,k_0);
    
    if systemproj
    % Use OSGeo/proj, when available. 
        
        prj4 = [ '-I +proj=tmerc +ellps=' opt.ellps ' +units' opt.units ];  
        prj4 = arrayfun(@(varargin) sprintf(...
            [prj4 ' +lat_0=%0.8f +lon_0=%0.8f +x_0=%0.3f +y_0=%0.3f +k_0=%0.6f'],varargin{:} ),...
            lat_0,lon_0,x_0,y_0,k_0,'unif',0);
        
        [Lon,Lat] = proj(prj4,x,y);
        return;     
    end

    if ~strcmpi(opt.ellps,'WGS84')
       warning('Custom ellipsoids require PROJ system package: using WGS84');
    end
    
    validateattributes(x,{'numeric'},{'real'},'proj','x');
    validateattributes(y,{'numeric'},{'real'},'proj','y');
    compatiblesize(x,y);
    
    % The following code was modified from: 
    % Rafael Palacios (2021). utm2deg, MATLAB Central File Exchange.
    % <https://www.mathworks.com/matlabcentral/fileexchange/10914-utm2deg>
    
    sa = 6378137.000000 ; sb = 6356752.314245;
    % e2 = (sa^2 - sb^2)/sb^2;
    % c = (sa^2) / sb;
    e2 = (sa/sb)^2 - 1;
    e = sqrt(e2);
    c = (sa/sb)*sa;
    alfa = 3/4 * e2;
    beta = 5/3 * alfa^2;
    gama = 35/27 * alfa^3;

    if any(lat_0 ~= 0)
        % [~,y_offset] = deg2utm(lat_0,lon_0,'lat_0',0,'lon_0',lon_0,'y_0',0);

        a1 = sind(2 * lat_0);
        a2 = a1 .* cosd(lat_0).^2;
        j2 = lat_0*pi/180 + (a1 / 2);
        j4 = (3 * j2 + a2) / 4;
        j6 = (5 * j4 + a2 .* cosd(lat_0).^2) / 3;
        y_offset = 0.9996 * c .* (lat_0*pi/180 - alfa .* j2 + beta .* j4 - gama .* j6);
        
        y_0 = y_0 - y_offset;
    end
    
    y = y - y_0;
    x = x - x_0;

    phi =  y / (6366197.724 * 0.9996);                                
    % v = c ./ sqrt(1 + e2 * cos(phi).^2) * 0.9996;
    % a = x ./ v;
    vi = hypot(1,e*cos(phi))./(0.9996*c);
    a = x .* vi;
    a1 = sin(2 * phi);
    a2 = a1 .* cos(phi).^2;
    j2 = phi + a1/2;
    j4 = (3*j2 + a2)/4;
    j6 = (5*j4 + a2.*cos(phi).^2)/3;
    
    Bm = 0.9996 * c .* (phi - alfa * j2 + beta * j4 - gama * j6);
    % b = (y - Bm) ./ v;
    b = (y - Bm) .* vi;
    Epsi = ((e2 * a.^2) / 2) .* cos(phi).^2;
    Eps = a .* (1 - (Epsi / 3));
    nab = b .* (1 - Epsi) + phi;
    sinheps = (exp(Eps) - exp(-Eps)) / 2;
    Delt = atan(sinheps ./ cos(nab));
    TaO = atan(cos(Delt) .* tan(nab));
    Lon = Delt *(180 / pi) + lon_0;
    Lat = phi + (1 + e2* cos(phi).^2 - 3/2*e2*sin(phi).*cos(phi).*(TaO - phi)).*(TaO - phi);
    Lat = Lat * (180 / pi);
end
