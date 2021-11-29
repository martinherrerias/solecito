function [lat,lon] = prj2abs(x,y,lat0,lon0,rotation,varargin)
% [LAT,LON] = PRJ2ABS(X,Y,LAT0,LON0,ROTATION) - Return latitude & longitude (WGS84 degrees)
%   given local project coordinates X, Y (meters).
%
%   The project coordinate system is a [rotated] project-centered transverse mercator projection,
%   whose origin lies at LAT0, LON0 (WGS84 degrees). Optional ROTATION is the angle (cw. degrees)
%   between geographic North and the Y+ axis.
%
%   If ROTATION is not provided, it will be assumed that projects in the Northern- (LAT0 > 0)
%   and Southern hemispheres use ROTATION 180° and 0°, respectively; so that their Y+ axis
%   points towards the Equator.
%
% See also: ABS2PRJ, UTM2DEG, DEG2UTM

    if nargin == 0, abs2prj(); return; end % common test

    % Assume convention -Y -> Equator
    if nargin < 5 || isempty(rotation), rotation = (lat0 < 0)*180; end
    assert(isnumeric(rotation) && isscalar(rotation) && isreal(rotation),'Bad rotation');

    assert(all(abs(lat0) <= 90) && all(abs(lon0) <= 180),...
        'Expecting decimal degree coordinates LAT0 and LON0');
    if all(abs(lat0) <= pi/2) && all(abs(lon0) <= pi)
        warning(['Unless your project is floating somewhere on the Gulf of Guinea, ',...
                 'consider checking your input LAT0, LON0 values']);
    end
    
    if rotation ~= 0
        [x,y] = compatiblesize(x,y);
        R = [cosd(rotation),-sind(rotation);sind(rotation),cosd(rotation)];
        XY = [x(:),y(:)]*R';
        x(:) = XY(:,1);
        y(:) = XY(:,2);
    end
 
    % Custom (project-centered) Transverse-Mercator Projection
    [lat,lon] = utm2deg(x,y,'lat_0',lat0,'lon_0',lon0,'x_0',0,'y_0',0,varargin{:});    
end
