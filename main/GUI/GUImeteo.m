function GUImeteo(varargin)
% GUIMETEO() - Imports Meteorological Data, calculates Sun Positions and derived quantities
%              Requires SimOptions, from GUIsetup
%
% GUIMETEO(MD,TIME,LOC) - Re-process existing METEODATA object.

    global GUI

    setflag('meteo',-3,'Running GUImeteo...');
    finalwillandtestament = onCleanup(@updateflags);
    fprintf('\nRunning GUImeteo...\n');
    
    opt = getSimOption('meteo');
    
    % Backdoor for re-filtering, resampling, etc. 
    if nargin > 0 && ~isa(varargin{1},'matlab.ui.control.UIControl')
        setflag('meteo',-3,'Reprocessing MeteoData...');
        fprintf('Reprocessing MeteoData\n');
        if nargin == 1 && isa(varargin{1},'MeteoData')
            MD = MeteoData(varargin{1});
        elseif nargin == 3 && all(cellfun(@isstruct,varargin))
        % GUImeteo(S,Time,Loc)
            MD = MeteoData(varargin{1},'t',varargin{2},'location',varargin{3},'interval','c');
        end
    else
        setflag('meteo',-3,'Importing Time & MeteoData...');
        fprintf('Importing Time & MeteoData\n');

        % Get Meteo Data
        filename = pickfile('*.meteo','Select a Meteo-Data file','fullpath',true);
        here = pwd(); lastwill = onCleanup(@() cd(here)); 
        cd(fileparts(filename)); % switch to meteo-file location for wps_.. files to be found
        
        fprintf('Source: %s\n',relativepath(filename));
        MD = MeteoData.import(filename,[]);
    end

    loctxt = @(S) [deg2dms(S.latitude,'NS'),', ',deg2dms(S.longitude,'EW')];

    if isfield(MD.location,'name')
        flagmsg{1} = sprintf('Site/dataset name: ''%s''',MD.location.name);
    else
        flagmsg{1} = 'Unnamed location';
    end
    if ischar(opt.filter), filtername = opt.filter; 
    elseif islogical(opt.filter), filtername = sprintf('logical, %d/%d samples',nnz(opt.filter),MD.Nt);
    else, filtername = sprintf('index, %d/%d samples',numel(opt.filter),MD.Nt);
    end
    assert(mod(opt.downsample,1)==0 && opt.downsample > 0,'Cannot use fractional re-sampling yet!');

    flagmsg{2} = sprintf('Loc: %s, %dmASL', loctxt(MD.location), MD.location.altitude);
    if ~isequal(filtername,'all')
        MD.info{end+1} = sprintf('filter:''%s''',filtername);
        flagmsg{2} = strjoin([flagmsg(2),MD.info(end)],', ');
    end
    if opt.downsample~=1
        MD.info{end+1} = sprintf('%dX downsampled',opt.downsample);
        flagmsg{2} = strjoin([flagmsg(2),MD.info(end)],', ');
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
                           'allow the calculation of POA irradiance.'],ds,loctxt(Loc),loctxt(Loc0));
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

    Nt0 = MD.Nt;

    setflag('meteo',-3,[flagmsg,{'Performing Quality Control & Ancilliary Gap Filling...'}]);
    MD = completemeteodata(MD);

    % PROVISIONAL
    newname = meteofilename(MD.location,MD.t(1),MD.t(end),MD.timestep,'suffix','_QCMD.mat');
    if isempty(dir(newname)), backupdelete(newname); end
    save(newname,'MD');
    
    MD = cleanup(MD);

    plot(MD,'shade');
    plot(MD,'ktrd');
    plot(MD,'heatmap');

    MD = MD.filter('downsample',opt.downsample,'filter',opt.filter);
    flagmsg{3} = sprintf('%d/%d time-steps, %s intervals\n',MD.Nt,Nt0,char(MD.timestep));
    %Nt = numel(MeteoData.GHI);
    fprintf('%s',flagmsg{3});

    % PROVISIONAL
    Loc = MD.location;
    Time = MD.t;
    [MD,SunPos] = legacy(MD);

    assignin('base','MD',MD);
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