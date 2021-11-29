function [R,ZXZ] = mountrotations(Mounts,sunaz,sunel,AzConv)
% R = MOUNTROTATIONS(MOUNTS,SUNAZ,SUNEL,[AZCONV]) - Returns a set of rotation matrices for the  
%   mount array represented by MOUNTS, at every solar position SUNAZ(t), SUNEL(t) with t = 1...Nt. 
% R = MOUNTROTATIONS(MOUNTS) - for 0a, returns the static rotation matrix for each tracker.
% [R,ZXZ] = MOUNTROTATIONS(...) - Additionally return equivalent intrinsic Euler-Angles ZXZ.
%
%   Rotation matrices R are to be left-multiplied by points [x';y';z'] defined in 'Mount-
%   Coordinate-System' (see note on coord. systems below). Points [x";y";0] on the active PV
%   surface map to points [X;Y;Z] in 'Project-Coordinate-System' using:
%
%       [x';y';z'] = [x";y";0] + Trackers.axisoffset;
%       [X;Y;Z] = R·[x';y';z'] + Trackers.centers
%
% INPUT:
%  SUNAZ : azimuth in convention Equator = 90°, positive around zenith. 
%       See SOLARPOSITION.FIXAZIMUTH to switch from other conventions.
%  SUNEL : sun elevation. 0° = sunset/sunrise, 90° = Zenith.
%
%  MOUNTS: structure defining a physical mount layout (see SAMPLESYSTEM for examples).
%   MOUNTS.type - char key {'2a','1aV','1aF','1aC','0a'} specifying tracking type: #a = #-axis,
%       V = vertical, F = Fix-tilt, C = Contour ("horizontal"). '0a' stands for static sheds.
%
%   MOUNTS fields {azimuth, tilt, slope, tracklimits, backtracking, groundcoverratio}
%       define mount/axis orientation and tracking behaviour, depending on MOUNTS.type:
%
%   - 2a: 2-axis trackers. Use only 'tracklimits' 4 vector [l r d u], where l, r are vertical-axis  
%            rotation limits and d, u are horizontal axis limits (0° = vertical, 90° = facing up).
%            Default is [-180 180 0 90] 
%
%   - 1aV: Vertica-axis trackers. Require common 'tilt' (scalar - plane angle vs ground) and use
%           'tracklimits' [l r], if available. Default is [-90,90]
%
%   - 1aC/1aF: One axis (tilted or horizontal) trackers. 
%        + tilt: (scalar/Nm-vector) - axis angle vs ground (default 0). More precisely, it is
%            the angle between the x'-(rotation)-axis, and the Project's Y-axis. Positive around X.
%        + azimuth: (scalar/Nm-vector) - horizontal angle between x'-axis and Project's Y-axis,
%            positive around Z. (default 0)
%        + tracklimits: rotation boundaries for roll around x' axis, default [-90,90]
%        + groundcoverratio: ratio between tracker width and row pitch (both perpendicular to
%            rotation axis), if backtracking == true, then groundcoverratio must be provided (*)
%
%       (*) PROVISIONAL: it is assumed that all trackers are controlled by an astronomical
%          algorithm, i.e. no tracker optimization to overcast conditions. All trackers rotate
%          by the same 'roll' angle, and ground-cover-ratio [GCR = sin(shade-free-angle)] is equal
%          for all trackers and east/west directions.
%          FUTURE: MOUNTS.group indices should allow the use of local, east-west GCR pairs.
%
%   - 0a: Fixed tables (sheds)
%        + tilt: (scalar/Nm-vector - approx. plane vs ground angle ('roll' around x' axis).
%        + azimuth: (scalar/Nm-vector, horizontal angle between x'-axis and Project's X-axis,
%            positive around Z. (default 0)
%        + slope: (Nt vector - vertical angle between x'- and X-axes.
%
% OUTPUT:
%  R: in the general case a 3·3·Nm*·Nt* array, where Rtrck(:,:,j,t) is the rotation matrix for 
%     tracker j at timestep t [SUNEL(t), SUNAZ(t)].
%     The array will be reduced if all matrices are equal along a given dimension:
%       [3,3,Nm,Nt] for trackers with different orientations
%       [3,3,1,Nt*] in all cases where all mounts share the same orientation at all times
%       [3,3,Nm*] for static mounts (0a, or others with equal rotation limits)
%       [3,3] for static, uniform mounts.
%
%  ZXZ: an Nt*·Nm*·3 array of equivalent intrinsic Euler-Angles in ZXZ convention. For mount j at
%     time t, rotation is equivalent to:
%           1. ZXZ(t,j,1) degrees around Z, followed by...
%           2. ZXZ(t,j,2) degrees around X', followed by...
%           3. ZXZ(t,j,3) degrees around Z"
%     See POLYGON3D.ROTMAT for further information.
%
% NOTE on coordinate systems: all coordinate systems obey the right hand rule: z = cross(x,y)
% 
%  TCS - Mount/Tracker-Coordinate-System [x',y',z']:
%   - The PV surface is parallel to the x'y' plane and has normal z'. Any tracker rotations take
%     place around O' (x'= y'= z'= 0), but the PV surface might offset by MOUNTS.axisoffset.
%   - In single-axis mounts (1aV,1aC,1aF), x' axis is parallel to the rotation axis, and points
%     MOUNTS.azimuth + MOUNTS.rotation degrees ccw. from North.
%     For 2a mounts, x' is the elevation-rotation axis (y' points up at sunrise/sunset).
%     For 0a mounts, y' is at the "upper" end of the tables (for tilt-angles > 0) and x' goes
%       in the direction of the table rows. For tilt = 0, y'points MOUNTS.azimuth + MOUNTS.rotation
%       degrees ccw. from North.
%   - Origins O' are located at MOUNTS.centers
%
%  PCS - Project-Coordinate-System [X,Y,Z] (see CHECKCOORDSYSTEM)
%   - Z = Zenith 
%   - Y = horizontal vector MOUNTS.azimuth + MOUNTS.rotation degrees CCW from North
%   - X = horizontal vector MOUNTS.azimuth + MOUNTS.rotation degrees CCW from East
%
% TODO: switch to quaternion representation
%
% See also: GUILAYOUT, PLOTTRACKERARRAY, SAMPLESYSTEM, CHECKCOORDSYSTEM
		
    if nargin == 0, test(); return; end

    parsestruct(Mounts,{'type'},'class','char','nonempty');
    [~,trcktype] = parselist(Mounts.type,{'2a','1aV','1aF','1aC','0a'},'mount type');
    
    parsestruct(Mounts,{},'opt',{'tilt','azimuth','slope'},'-n','-r','-v','-c','size',[NaN,1]);
    if trcktype < 5
        parsestruct(Mounts,{},'opt',{'tracklimits'},'-n','-r','-v',@(x) any(numel(x) == [2,4]));
        parsestruct(Mounts,{},'opt',{'backtracking'},'class',{'numeric','logical'},'binary');
        parsestruct(Mounts,{},'opt',{'groundcoverratio'},'-n','-r','-s');
    end
    
    % Check track limits
    if trcktype ~= 5
        if ~isfield(Mounts,'tracklimits'), trcklims = Mounts.tracklimits(:)';
        else, trcklims = [-180 180 0 90];
        end
        if ~isfield(Mounts,'tracklimits') || (trcktype == 1 && numel(trcklims) < 4)
            if (trcktype == 1 && numel(trcklims) < 4), trcklims(3:4) = [0,90]; end
            warning('mountrotations:nolims','Rotation limits not defined: Using no limits!');
        end
        assert(trcklims(2) >= trcklims(1) && trcklims(4) >= trcklims(3),...
            'Expecting tracking limits [Hm,HM,Vm,VM] with HM >= Hm [and Vm >= VM]');
        
        % Find out if the system can actually move i.e. ~all(constrainedaxis)
        constrainedaxis = [trcklims(1)==trcklims(2),trcktype ~= 1 || (trcklims(3)==trcklims(4))];        
    else
       constrainedaxis = [true,true];
    end

    % Check solar elevation and azimuth angles
    if nargin < 3
        if all(constrainedaxis), sunaz = 0; sunel = 0; % not like it matters
        else
            error('mountrotations:nargs',...
                'Solar azimuth and elevation required for calculation.');
        end
    else
        validateattributes(sunaz,{'numeric'},{'real','vector'});
        validateattributes(sunel,{'numeric'},{'real','vector','size',size(sunaz)});
        if size(sunaz,2) > 1, sunaz = sunaz'; sunel = sunel'; end
        
        if nargin < 4 || isempty(AzConv), AzConv = 'N2E'; end
        if ~isequal(AzConv,'N2E')
            sunaz = solarposition.fixazimuth(sunaz,AzConv,'N2E',Mounts.origin(2));
        end
    end

    % all trackers except 2a and 1aV require azimuth
    if trcktype > 2 
        if ~isfield(Mounts,'azimuth') || isempty(Mounts.azimuth)
            azimuth = 0; 
            warning('mountrotations:NoAzimuth','Azimuth zero (1a-EW/0a-equator) will be assumed')
        else
            azimuth = Mounts.azimuth;
        end
    end
    
    % Adjust solar azimuth to project convention: X = 0, Y = 90
    if ~isfield(Mounts,'rotation'), Mounts.rotation = 0; end
    sunaz = 90 - sunaz - Mounts.rotation;
    
    % all trackers except 2a require tilt
    if trcktype > 1 
        if ~isfield(Mounts,'tilt') || isempty(Mounts.tilt)
            tilt = 0; 
            warning('mountrotations:NoTilt','Tilt zero (horizontal) will be assumed')
        else
            tilt = Mounts.tilt; 
        end
    end
    
    % backtracking parameters for 1aC and 1aF
    if any(trcktype == [3,4]) && ~constrainedaxis(1) 
        if ~isfield(Mounts,'backtracking') || isempty(Mounts.backtracking)
            backtracking = false;
            warning('mountrotations:Backtracking','No-backtracking assumed!')
        else
            backtracking = Mounts.backtracking;
        end
        if backtracking 
            if ~isfield(Mounts,'groundcoverratio') || isempty(Mounts.groundcoverratio)
                error('mountrotations:backtracking requires knowledge of ground-cover-ratio!');
            else
                gcr = Mounts.groundcoverratio; 
                assert(gcr >= 0 && gcr <= 1,'mountrotations:gcr','un-physical ground-cover-ratio');
            end
        else
            gcr = 0;
        end
    else
       backtracking = false; gcr = 0; % just for parsing, won't be used
    end

    % Do the real work
    switch Mounts.type
        case '0a'
            [R,ZXZ] = rotations_0a(tilt,azimuth,Mounts.slope);
                
        case '2a'
            [R,ZXZ] = rotations_2a(sunaz,sunel,trcklims);
            
        case '1aV' % simulate as constrained 2a system (*)
            trcklims(3:4) = 90 - Mounts.tilt;
            [R,ZXZ] = rotations_2a(sunaz,sunel,trcklims);
            return
  
        case {'1aF','1aC'}
            [R,ZXZ] = rotations_1aX(sunaz,sunel,trcklims,azimuth,tilt,backtracking,gcr);
    end
    
    %if nargout > 2
    % Get an approximate mean-rotation matrix for each timestep
    % Ravg = squeeze(mean(R,3));
    % zxz0(:,1) = atan2d(Ravg(1,3,:),-Ravg(2,3,:));
    % zxz0(:,2) = atan2d(hypot(Ravg(1,3,:),Ravg(2,3,:)),Ravg(3,3,:));
    % zxz0(:,3) = atan2d(Ravg(3,1,:),Ravg(3,2,:));
    % Ravg = polygon3d.rotmat(zxz0,'ZXZ');
end

function [R,ZXZ] = rotations_0a(tilt,azimuth,slope)
% Calculates rotation matrixes for an array of fixed mounts (0a trackers) originally laying flat
% on the ground with x' = X, y' = Y

    [tilt,azimuth,slope] = compatiblesize(tilt,azimuth,slope);
    Nm = max([numel(tilt),numel(azimuth),numel(slope)]);

    equaltilt = size(unique([tilt(:),azimuth(:),slope(:)],'rows'),1) == 1;
    if Nm > 1 && equaltilt      
    % reduce if all share the same orientation
        tilt = tilt(1); azimuth = azimuth(1); slope = slope(1);
    end

    [R,ZXZ] = polygon3d.rotmat([azimuth(:),-slope(:),tilt(:)],'ZYX');
    ZXZ = permute(ZXZ,[3,1,2]);
end

function [R,ZXZ] = rotations_2a(sunaz,sunel,trcklims)
% Rotation matrixes for a 2a tracker originally laying flat on the ground with x' = X, y' = Y.
% sunaz must be the CCW angle from X! i.e. E2N rotated by MOUNTS.rotation

    phi = sunaz(:) + 90;
    phi = solarposition.fixtoplusminus180(phi);
    phi = min(max(phi,trcklims(1)),trcklims(2));  % bind to vertical axis limits

    theta = 90-sunel(:);                                    % rotation around x'...
	theta  = max(min(theta,90-trcklims(3)),90-trcklims(4)); % bound to limits
	
    ZXZ = [phi,theta,zeros(size(phi))];     % [Nt,3]
    R = polygon3d.rotmat(ZXZ,'ZXZ');        % [3,3,Nt]
    R = permute(R,[1 2 4 3]);               % [3,3,1,Nt]
    ZXZ = permute(ZXZ,[1,3,2]);             % [Nt,1,3]
end

function [R,ZXZ] = rotations_1aX(sunaz,sunel,trcklims,azimuth,tilt,backtracking,gcr)
% Calculates rotation matrices for an array of 1aC trackers, originally laying flat on the ground
%   with with x' = X, y' = Y.
% sunaz must be the CCW angle from X! i.e. E2N rotated by MOUNTS.rotation
 
    Nm = max(numel(tilt),numel(azimuth));
    if Nm>1 && isscalar(tilt),tilt(1:Nm) = tilt; end
    if Nm>1 && isscalar(azimuth),azimuth(1:Nm) = azimuth; end
    
    equaltilt = isscalar(tilt) || all(tilt==tilt(1));
    equalazimuth = isscalar(azimuth) || all(azimuth == azimuth(1));
    
    if Nm > 1 && equaltilt && equalazimuth      
    % reduce if all share the same orientation
        Nm = 1; tilt = tilt(1); azimuth = azimuth(1);
    else
        tilt = tilt(:)'; azimuth = azimuth(:)';
    end
    
    % Get preliminary rotation matrix (to get average x' and z')
    Rtrck0 = polygon3d.rotmat([mean(azimuth)+90,-mean(tilt)],'ZY');
        
    fixed = trcklims(1)==trcklims(2);   % Stuck (disguised 0a case)
    if fixed
        Nt = 1;                 % reduce time-dimension
        roll = trcklims(1);     % constant roll-angle
    else
        Nt = numel(sunel);
        if ~equalazimuth
            warning('mountrotations:AzVector','Average azimuth will be used for tracking')
        end

        % list of unitary vectors with the direction of solar beams (downwards)
        sunvec = [-cosd(sunaz').*cosd(sunel');-sind(sunaz').*cosd(sunel');-sind(sunel')];

        % roll = atan(s·y'/-s·z'), an Nt vector
        roll = atan2(sunvec'*Rtrck0(:,2),-sunvec'*Rtrck0(:,3))*180/pi;

        % Correct for backtracking: r = r* + sign(r*)·(asin(cos(r*)·a/b)-90°)   
        % Note: comes from Law of Sines: sin(90°-r*+r)/a = sin(90°-r*)/b
        %       nicer looking r = r* - sign(r*)·acosd(cos(r*)·a/b) doesn't work
        if backtracking
            roll = min(max(roll,-90),90);
            iout = cosd(roll) < gcr; % sun below shade-free angle
            roll(iout) = sign(roll(iout)).*(asind(cosd(roll(iout))/gcr)-90)+roll(iout);
            roll(sunel < 0) = 0;
        end

        % Clip to tracking limits
        roll = min(max(roll,trcklims(1)),trcklims(2));
    end
    
    azimuth = repmat(azimuth,Nt,1);
    tilt = repmat(tilt,Nt,1);
	
	R = zeros(3,3,Nm,Nt);
    ZXZ = zeros(Nt,Nm,3);
    for tr = 1:Nm
        [r,zxz] = polygon3d.rotmat([azimuth(:,tr)+90,-tilt(:,tr),roll],'ZYX');
        R(:,:,tr,:) = permute(r,[1 2 4 3]);
        ZXZ(:,tr,:) = permute(zxz,[1,3,2]);
    end
end

function test()

    h = GUIfigure('mountrotations','samplesystem/mountrotations test','2:1');
    ax = subplot(1,3,1:2);
    ax(2) = subplot(1,3,3);
    ax(2).Visible = false;
    p = uipanel(h,'Position',ax(2).Position);
    
    btn = uicontrol('Parent',p,'Style','pushbutton','String','Animate',...
        'Units','normalized','Visible','on','Position',[0.1 0.1 0.8 0.1]);
    
    H = plotcontrols({'m','s','s','s'},{'Mount Type','Latitude','Azimuth','Rotation'},...
        {{'0a','1aF','1aC','2a','1aV'},[-80,80,1,5],[-90,90,1,5],[-180,180,1,5]},...
        {'0a',25,0,15},@update,'parent',p,'position',[0.1 0.4 0.8 0.6],'-skipupdate');
    set(H,'BusyAction','cancel');

    update('0a',25,0,15);
    view([300,20]);

    function update(type,lat,azimuth,rotation,~,evt)
        btn.UserData.running = false; % stop animation
        if nargin > 4 && ~strcmp(evt.EventName,'Action'), return; end
        
        set(H,'Enable','off');
        % p.Enable = 'off';
        Loc = struct('latitude',lat,'longitude',0,'altitude',0,'TimeZone',0);
        [az,el] = view();
        Trck = samplesystem(type,[],'location',Loc,'landscape',true,'az',azimuth,...
            'rotation',rotation);
        plottrackerarray(Trck,'btn',btn,'view',[az,el],ax(1));
        % view(gca,0,90);
        title(sprintf('%s, lat = %d°',type,lat));
        % p.Enable = 'on';
        set(H,'Enable','on');
    end
end
