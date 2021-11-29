function GUIterrain(~,~)
% Read DEM file, if available, and calculate horizon profiles
% Requires SimOptions.flathorizon, from GUIsetup
% Requires Location.latitude if ~flathorizon, from GUImeteo

    setflag('terrain',-3,'Runing GUIterrain...');
    finalwillandtestament = onCleanup(@updateflags);
    
    fprintf('\nRunning GUIterrain...\n');
    flathorizon = getSimOption('flathorizon');
    flatthreshold = getSimOption('terrain.flatthreshold');
    
    Trackers = evalin('base','Trackers');

    printstr = @(s) fprintf('\t%s\n',s);

    HP = struct('nHor',[],'fHor',[]); % if nothing happens, this will be it
    if flathorizon          
    % Issue warnings if flathorizon flag is set, skip the rest
        if ~isempty(pickfile({'*.tif','*.xyz','*.hor'},Inf)) % terrain data is available 
            setflag('terrain',1,{' Flat horizon! - type SimOptions.flathorizon = false',...
                                 ' at console to use available *.hor/*.xyz/*.tif file(s)'});
        else
            setflag('terrain',1,{' Flat horizon - type SimOptions.flathorizon = false',...
                                 ' at console to browse for *.hor/*.xyz/*.tif file(s)'});
        end

        ang = atan2d(Trackers.centers(3,:),hypot(Trackers.centers(1,:),Trackers.centers(2,:)));
        if max(abs(ang)) > flatthreshold || ~isempty(pickfile({'*.tif','*.xyz','*.hor'},Inf))
            warning('GUIterrain:notflat',['Using Flat-Horizon option on possibly complex '...
                'terrain, use SimOptions.flathorizon = false for a realistic calculation']);
        else
            printstr('Using Flat Horizon for all mounts');
        end
        assignin('base','HorizonProfile',HP);
        return
    end

    hpfilename = ['horprof_' Trackers.name '.mat'];
    if isempty(dir(hpfilename))
        hpfilename = fullfile(fileparts(getSimOption('prjname')),hpfilename);
    end

    if ~isempty(dir(hpfilename))
        switch optquestdlg('Near-Horizons seem already available, do you want to use them?')
        case 'Yes'
            load(hpfilename,'HP');
        % checklocation(Loc);
            if numel(HP.nHor) ~= size(Trackers.centers,2)
                flagmsg{1} = 'Pre-calculated results could not be used';
                HP = struct('nHor',[],'fHor',[]);
            else
                flagmsg{1} = 'Using precalculated Near-Horizon-Profiles';
            end
        case 'No'
            flagmsg{1} = 'Not using precalculated Near-Horizon-Profiles';
        otherwise
            error('GUIterrain:questdlg','Cancelled by user')
        end
    else
        flagmsg{1} = 'No precalculated terrain results found';
    end
    printstr(flagmsg{1});

    if isempty(HP.nHor)
    % Read DEM from geotiff/xyz file(s), and calculate near and far horizons
    
        % Get circumcenter and circumradius of the plant
        [R0,C] = MinBoundCircle(Trackers.centers(1:2,:)');
        % Get the WGS84 coordinates of C
        [lat_C,lon_C] = checkcoordsystem(C(1),C(2),0,Trackers,'output','abs');

        % Two possible criteria to define the minimum required span of DEM data:
        R1 = getSimOption('terrain.demradius');  % 1. Min. radius for far-horizon determination
        R2 = R0/tand(getSimOption('angtol'));    % 2. Min. radius for meaningful near-horizons

        % If a far-horizon is available, the second is important, otherwise the first
        if ~isempty(dir('*.hor')), R0 = min(R2,R1); else, R0 = max(R2,R1); end
        setSimOption('terrain.demradius',R0);

        demfiles = {};
        switch numel(dir('*.tif'))
        case 0
        % No TIFF files found...
            switch optquestdlg(['No DEM-raster (*.tif) found, do you want to browse for a ',...
                'file somewhere else, attempt to generate one from ASTER V2 tiles, or ',...
                'ignore far terrain and get local-horizons?'],'Missing DEM raster',...
                'Browse','Generate','Ignore','Browse')

                case 'Browse'
                % Go fetch it somewhere else
                    demfiles = {pickfile('*.tif')};

                case 'Generate'
                % If no DEM raster is available, attempt to download and merge from ASTER
                    demfiles = {getgeotiff(lat_C,lon_C,R0)}; 

                case 'Ignore' % go on as if nothing happened
                otherwise, error('Stopped by user');
            end
        case 1
        % Single file found, try to use it directly
            demfiles = pickfile('*.tif',Inf); % cell-array with string(s)
        otherwise
        % If no DEM raster is available, attempt to download and merge from ASTER
            demfiles = {getgeotiff(lat_C,lon_C,R0)}; 
        end
        if ~isempty(dir('*.xyz')), demfiles{end+1} = pickfile('*.xyz'); end

        if ~isempty(demfiles)
            setflag('terrain',-3,'Reading geotif file...');
            Terr = readterrain(Trackers,demfiles);
            flagmsg{end+1} = sprintf('%0.2e - point terrain model imported successfully',numel(Terr.x));
            setflag('terrain',-3,flagmsg); printstr(flagmsg{end});

            % Calculate near- and far-horizons from DEM
            setflag('terrain',-3,[flagmsg,{'Calculating horizons from DEM...'}]);
            HP = gethorizons(Trackers,Terr);
            flagmsg{end} = sprintf('Horizons calculated from %0.2e-point DEM',numel(Terr.x));
            setflag('terrain',-3,flagmsg); printstr(flagmsg{end});

            save(hpfilename,'HP');       
        else
        % If no geotiff file is available, do your best without it (use approxgroundmesh internally)
            flagmsg{1} = 'No DEM data found, immediate-horizons from mount coordintes'; 
            HP = gethorizons(Trackers);
            setflag('terrain',-3,flagmsg); printstr(flagmsg{1});
        end
    end

    if ~isempty(dir('*.hor'))
    % With or without DEM, read far-horizon directly from *.hor file
        flagmsg{end+1} = 'Reading and merging Far-Horizon Profile from *.hor file';
        setflag('terrain',-3,flagmsg); printstr(flagmsg{end});

        [az,el] = readhorizon(pickfile('*.hor'));
        az = solarposition.fixazimuth(az,'S2W','E2N') - Trackers.rotation;  % project convention	
        hor = polygon3d(az,el,1,'sph');

        fprintf('\tImported Far-Horizon: %d points, %0.1f° mean elevation\n',...
                            numel(az),(el'*az([2:end,1])-az'*el([2:end,1]))/720);
                        
        % Merge with previously calculated horizons
        prj = polyprojector('stereo');
        hor = prj.project(hor,[0;0;1],true);
        q = prj.project(HP.fHor,[0;0;1],true);
        q = intersectpolygons(hor,q);
        HP.fHor = prj.inverse(q);
        [az,el] = cart2sph(HP.fHor.x,HP.fHor.y,HP.fHor.z);
        
        for j = 1:numel(HP.nHor)
            q = prj.project(HP.nHor(j),[0;0;1],true);
            q = intersectpolygons(hor,q);
            HP.nHor(j) = prj.inverse(q);
        end

        flagmsg{end} = sprintf('Merged Far-Horizon: %d points, %0.1f° mean elevation',...
                        numel(az),(el*az([2:end,1])'-az*el([2:end,1])')/(4*pi)*180/pi);
        setflag('terrain',-3,flagmsg); 
        printstr(flagmsg{end});

        % Plot in existing figure from gethorizons
        GUIfigure('horizons');
        plot3(HP.fHor.x([1:end,1]),HP.fHor.y([1:end,1]),HP.fHor.z([1:end,1]),'k-');
    end

    flathorizon = all([HP.fHor.z] == 0);
    for j = 1:numel(HP.nHor)
        if ~flathorizon, break; end
        flathorizon = max(HP.nHor(j).z) < flatthreshold;
    end
    if flathorizon
        setSimOption('flathorizon',true);
        flagmsg{end+1} = 'Everything looks flat - setting opt.flathorizon';
        setflag('terrain',-3,flagmsg); 
        printstr(flagmsg{end});
        printstr('Use SimOptions.terrain.flatthreshold = 0 to avoid this simplification');
        setflag('horizon',1,' Flat horizon read from *.hor file');
    end

    assignin('base','HorizonProfile',HP);

    % If it went without errors, move on
    setflag('terrain',1,flagmsg);
end

function [az,el] = readhorizon(filename)
    fileID = fopen(filename);
    if fileID < 0, error('readsplitarraydef:fopen','Could not open file'); end
    firstlines = textscan(fileID,'%[^\n\r]',2);
    fclose(fileID);

    % Guess delimiters from second line
    delim = [char(9),';, ']; % delimiter priority: (tab) ; , (space)
    for j = 1:numel(delim)
        if any(firstlines{1}{2}==delim(j))
            delim = delim(j);
            break
        end
    end

    % Try to split and read a number
    %fields = strsplit(firstlines{1}{1},delim);
    %skipheader = isnan(str2double(fields{1}));

    hordata = dlmread(filename,delim,1,0);
    hordata = sortrows(hordata,1);
    az = hordata(:,1); el = hordata(:,2);
end

function checklocation(Loc)
% Check acquired location against that of Meteo-Data, if available

    global GUI;
    if ~(GUI.meteo.flag > 0), return; end
    
    Loc0 = evalin('base','Location'); 
    ds = solarposition.arcdist(Loc.latitude,Loc.longitude,Loc0.latitude,Loc0.longitude);
    ds = ds * 6370; % aprox distance in km
    if ds > getSimOption('maxlocdist')
        msg = sprintf(['The coordinates for the Meteo-Data you are using seem '...
                       '%0.1f km away from the project. The coordinates used for '...
                       'terrain calculations (%s) will not match those used for '...
                       'irradiance-transposition (%s), which will likely not ' ...
                       'allow the calculation of POA irradiance.'],ds,txt(Loc),txt(Loc0));
        switch optquestdlg([msg,'Do you still want to continue?'],...
                'Inconsistent coordinates',...
                'No, let me fix it','Yes, see what happens','No, let me fix it')
        case 'Yes, see what happens'
            warning([msg '!']);
        otherwise
            setflag('terrain',0,{'Stopped to fix to inconsistent coordinates',...
                                 'Ready when you are'});
            return;
        end
    end
end
