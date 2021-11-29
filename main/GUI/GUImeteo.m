function GUImeteo(varargin)
% GUIMETEO() - Imports Meteorological Data, calculates Sun Positions and derived quantities
%              Requires SimOptions, from GUIsetup
%
% GUIMETEO(METEODATA,TIME,LOC) - Re-process existing METEODATA.

    global GUI

    setflag('meteo',-3,'Running GUImeteo...');
    finalwillandtestament = onCleanup(@updateflags);
    fprintf('\nRunning GUImeteo...\n');
    
    opt = getSimOption('meteo');
    
    % Backdoor for re-filtering, resampling, etc. 
    if nargin == 3 && all(cellfun(@isstruct,varargin))
        setflag('meteo',-3,'Reprocessing MeteoData...');
        fprintf('Reprocessing MeteoData\n');
        [MeteoData, Time, Loc] = deal(varargin{:});
        [Time,MeteoData.timestep] = parsetime(Time);
    else
        setflag('meteo',-3,'Importing Time & MeteoData...');
        fprintf('Importing Time & MeteoData\n');

        % Get Meteo Data
        filename = pickfile('*.meteo','Select a Meteo-Data file','fullpath',true);
        here = pwd(); lastwill = onCleanup(@() cd(here)); 
        cd(fileparts(filename)); % switch to meteo-file location for wps_.. files to be found
        
        fprintf('Source: %s\n',relativepath(filename));
        [MeteoData, Time, Loc] = getMeteoData(filename,[]);
    end
    Nt = numel(Time);

    txt = @(S) [deg2dms(S.latitude,'NS'),', ',deg2dms(S.longitude,'EW')];

    if isfield(Loc,'name')
        flagmsg{1} = sprintf('Site/dataset name: ''%s''',Loc.name);
    else
        flagmsg{1} = 'Unnamed location';
    end
    if ischar(opt.filter), filtername = opt.filter; 
    elseif islogical(opt.filter), filtername = sprintf('logical, %d/%d samples',nnz(opt.filter),Nt);
    else, filtername = sprintf('index, %d/%d samples',numel(opt.filter),Nt);
    end
    assert(mod(opt.downsample,1)==0 && opt.downsample > 0,'Cannot use fractional re-sampling yet!');

    flagmsg{2} = sprintf('Loc: %s, %dmASL', txt(Loc), Loc.altitude);
    if ~isequal(filtername,'all')
        MeteoData.info{end+1} = sprintf('filter:''%s''',filtername);
        flagmsg{2} = strjoin([flagmsg(2),MeteoData.info(end)],', ');
    end
    if opt.downsample~=1
        MeteoData.info{end+1} = sprintf('%dX downsampled',opt.downsample);
        flagmsg{2} = strjoin([flagmsg(2),MeteoData.info(end)],', ');
    end
    fprintf('%s\n',flagmsg{:});

    if GUI.terrain.flag > 0 && evalin('base','isfield(Trackers,''Location'')')
    % Check acquired location against that of GUIterrain, if available
        Loc0 = evalin('base','Trackers.Location'); 
        ds = solarposition.arcdist(Loc.latitude,Loc.longitude,Loc0.latitude,Loc0.longitude);
        ds = ds * 6370; % aprox distance in km
        if ds > getSimOption('maxlocdist')
            msg = sprintf(['The coordinates for the Meteo-Data you are using seem '...
                           '%0.1f km away from the project. The coordinates used for '...
                           'irradiance-transposition (%s) will not match those used for '...
                           'shading-calculations (%s), which might not ' ...
                           'allow the calculation of POA irradiance.'],ds,txt(Loc),txt(Loc0));
            switch optquestdlg([msg,'Do you still want to continue?'],...
                    'Inconsistent coordinates',...
                    'No, let me fix it','Yes, see what happens','No, let me fix it')
            case 'Yes, see what happens'
                warning([msg '!']);
            otherwise
                setflag('meteo',0,{'Stopped to fix to inconsistent coordinates',...
                                     'Ready when you are'});
                return;
            end
        end
    end
    
    % TODO: Get sensor data, cluster variables by type
    % [MeteoData,Sensors] = getsensordata(MeteoData);
    % ...

    % Calculate Solar-Position & dependent parameters
    setflag('meteo',-3,[flagmsg,{'Calculating solar position & dependent variables...'}]);
    [MeteoData, Time, SunPos, Loc] = completemeteodata(MeteoData,Time,Loc);
    flagmsg{3} = sprintf('%d/%d time-steps, %s intervals\n',...
        numel(Time),Nt,char(MeteoData.timestep));
    %Nt = numel(MeteoData.GHI);
    fprintf('%s',flagmsg{3});

    % struct2csv({Time,SunPos,MeteoData},Nt,'meteo.csv');

    assignin('base','MeteoData',MeteoData);
    assignin('base','Time',Time);
    assignin('base','SunPos',SunPos);
    assignin('base','Location',Loc);
    
    % If it went without errors, move on
    setflag('meteo',1,flagmsg);

    % If Array-Definition and Electrical-Models also exist, check them with new MeteoData
    if GUI.arrdef.flag > 0 && GUI.models.flag > 0
       ok = checkdesign('-verbose','-soft');
       if ~ok, setflag('models',2); else, setflag('models',1); end
    end
end
