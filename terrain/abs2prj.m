function [x,y] = abs2prj(lat,lon,lat0,lon0,rotation,varargin)
% [X,Y] = ABS2PRJ(LAT,LON,LAT0,LON0,ROTATION) - Return project coordinates X, Y (meters), given
%   latitude & longitude (WGS84 degrees)
%
%   The project coordinate system is a [rotated] project-centered transverse mercator projection,
%   whose origin lies at LAT0, LON0 (WGS84 degrees). Optional ROTATION is the angle (cw. degrees)
%   between geographic North and the Y+ axis.
%
%   If ROTATION is not provided, it will be assumed that projects in the Northern- (LAT0 > 0)
%   and Southern hemispheres use ROTATION 180° and 0°, respectively; so that their Y+ axis
%   points towards the Equator.
%
% See also: PRJ2ABS, UTM2DEG, DEG2UTM

    if nargin == 0, test(); return; end
    narginchk(4,6);

    % Assume convention -Y -> Equator (so that for 0a, positive tilts face equinox noon)
    if nargin < 5 || isempty(rotation), rotation = (lat0 < 0)*180; end
    
    [lat,lon] = compatiblesize(lat,lon);
    assert(~any(abs([lat(:);lat0]) > 90) && ~any(abs([lon(:);lon0]) > 180),...
        'Expecting decimal degree coordinates LAT, LAT0, LON, LON0');
    if all(abs([lat(:);lat0]) <= pi/2) && all(abs([lon(:);lon0]) <= pi)
        warning(['Unless your project is floating somewhere on the Gulf of Guinea, ',...
                 'consider checking your input LAT, LAT0, LON, LON0 values']);
    end
    
    % Custom (project-centered) Transverse-Mercator Projection
    [x,y] = deg2utm(lat,lon,'lat_0',lat0,'lon_0',lon0,'x_0',0,'y_0',0,varargin{:});
    
    if rotation ~= 0
        R = [cosd(rotation),-sind(rotation);sind(rotation),cosd(rotation)];
        XY = [x(:),y(:)]*R;
        x(:) = XY(:,1);
        y(:) = XY(:,2);
    end

end

function test()
    x = rand(10,1)*1000; 
    y = rand(10,1)*1000;
    
    r = rand(1)*360-180;
    % r = [];
    lat0 = round(rand(1)*160-80,2); 
    lon0 = round(rand(1)*360-180,2);
    
    [lat,lon] = prj2abs(x,y,lat0,lon0,r);
    [x2,y2] = abs2prj(lat,lon,lat0,lon0,r);
    
    ax = [1000 0 0; 0 0 1000]';
    
    GUIfigure('abs2prj.test','prj2abs + abs2prj round test','2:1'); clf(); 
    h = arrayfun(@(j) subplot(1,2,j),1:2);
    for j = 1:numel(h)
        axis(h(j),'equal');
        grid(h(j),'on');
        hold(h(j),'on');
    end
    
    title(h(1), sprintf('%0.2f°N, %0.2f°E, r = %0.1f°',lat0,lon0,r))
    plot(h(1),x,y,'ro')
    plot(h(1),x2,y2,'bx');
    plot(h(1),ax(:,1),ax(:,2),'c');
    xlabel(h(1),'x (m)')
    ylabel(h(1),'y (m)')
    
    plot(h(2),lon,lat,'ro');
    xlabel(h(2),'lon')
    ylabel(h(2),'lat')
    
    [ax(:,2),ax(:,1)] = prj2abs(ax(:,1),ax(:,2),lat0,lon0,r);
    plot(h(2),ax(:,1),ax(:,2),'c');
end