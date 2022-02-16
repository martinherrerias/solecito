function GUIshading(~,~)
% GUISHADING() - Search for precalculated shading-results, or perform new shading analysis. 
%
% For existing results, recognized file-names will have the forms: 
%
%   PATH/LAYOUTNAME_TYPE.shdres or
%   PATH/PRJNAME_TYPE.shdres
%
%   Where PATH can be either PWD or the project's path, and TYPE can take the forms:
%
%   xGy = Full year, grid-based analysis, x.y° resolution (e.g. 1G25 for 1.25°)
%   xIy = Partial/irregular grid-based analysis, x.y° resolution
%   nP = n particular points step-by-step analysis
%
% If several candidate-files are found, priority is G > I > P, in ascending resolution.
%
% For new analysis, grid-based analysis will be preferred whenever:
%
%   - Meteo-Data (and solar-positions) are not yet available
%   - The number of simulation steps is higher than the required grid steps, or...
%   - The user has decided so (upon UI dialog)
%
%   In any of these cases, the grid resolution is controlled by SimOptions.shadingresolution, and
%   the result will be stored on PRJPATH/LAYOUTNAME_xGy.shdres or PRJPATH/LAYOUTNAME_xIy.shdres,
%   depending on whether a full-year (known location) or partial-grid (from solar-positions) was
%   used.
%
%   Otherwise, step-by-step analysis will be saved as PRJPATH/PRJNAME_nP.shdres

    global SimOptions
    global GUI
    
    SimOptions.groundshading = SimOptions.groundshading && SimOptions.diffuseshading;
    SimOptions.anisotropic = SimOptions.anisotropic && SimOptions.diffuseshading;
   
    setflag('shading',-3,'Runing GUIshading...');
    finalwillandtestament = onCleanup(@updateflags);
    
    fprintf('\n\nRunning GUIshading...\n');
        
    % delete(findall(0,'Tag','TMWWaitbar'));
    % h = GUIfigure('plotter');
    % if ishandle(h), delete(h); end % otherwise MATLAB crashes

    Trck = evalin('base','Trackers');
    HorProf = evalin('base','HorizonProfile');
    
    if ~isfield(Trck,'analysedtrackers'), Trck.analysedtrackers = 1:size(Trck.centers,2); end
    if ~isfield(Trck,'masks')
    % Calculate tracker masks
        masktolerance = SimOptions.masktolerance;
        fprintf('\nCalculating tracker masks (%0.1e mask tolerance)...\n',masktolerance);
        Trck.masks = gettrackermasks(Trck,masktolerance);
        fprintf('%0.1f (avg.) neighbors/tracker\n',mean(sum(Trck.masks,2)));
    end

    % Get ShadingRegions (sky-tessellation) from options
    SR = ShadingRegions('auto');
    SR = rotate(SR,-Trck.rotation);
    flagmsg{1} = ['Geometrical framework: ',SR.description];
    
    % Get IAM from base/options
    [~,flagmsg{2},IAM] = checkIAM();
    fprintf('\t%s\n',flagmsg{:});
    
    % Look for existing Shading-Results (project | layout-specific)
    [ShRes,flagmsg{3}] = lookforprevious(Trck,SR);
    frombackup = ~isempty(ShRes) && ~isa(ShRes,'ShadingResults');
   
    if ~isempty(ShRes) && ~frombackup
        assignin('base','ShRes',ShRes);
        setflag('shading',1,flagmsg);
        return
    end
    setflag('shading',-3,flagmsg);
  
    % Run/resume shading analysis  
    
    if frombackup
        resfilename = ShRes.options.backup;
    else
        % Decide on whether to do full-year(interpolation) or a partial (point-wise) analysis
        [SunPos,type,flagmsg{4}] = analysistype(SimOptions);
        if contains(type,'P')
            resfilename = shadingfilename(type); % SimOptions.prjname
        else
            resfilename = shadingfilename(type,Trck.name);
        end
        fprintf('\nRunning Shading Analysis...\n');
        setflag('shading',-3,[flagmsg,{'Runing Shading Analysis...'}]);
    end

    parallel = SimOptions.runparallel && ~frombackup;
    % if isequal(Trck.type,'0a')
    %     opt.runparallel = false;
    % end
    if parallel
        M = getnworkers([],'vars',whos('SunPos','Trck','HorProf','SR'));
        if M <= 1
            warning('Not enough memmory to span multiple workers');
            parallel = false; 
        end
    end
    
    if frombackup
        ShRes = ShadingAnalysis(SunPos,Trck,HorProf,SR,'IAM',IAM,'backup',ShRes);
    elseif parallel
        % Run shading analysis on parallel cluster
        [ShRes,cleaner] = runparallel(@ShadingAnalysis,{SunPos,Trck,HorProf,SR,'IAM',IAM},1,...
            [],'N',numel(SunPos.El),'-simoptions','-backup');
    else
        % DEBUG Backdoor: Run single-threaded
        ShRes = ShadingAnalysis(SunPos,Trck,HorProf,SR,'IAM',IAM,'backup',[resfilename,'~']);
    end
    
    % Add info fields
    ShRes.info.type = type;
    ShRes.info.layout = Trck.name;
    
    save(resfilename,'ShRes');

    timestr =@(x) [num2str(floor(x/3600)) ':' datestr(x/(24*3600),'MM:SS')];
    fprintf('\nShading analysis complete (%s elapsed)\n',timestr(ShRes.info.timing.globaltime));
    flagmsg{end+1} = 'Shading analysis completed successfully';      

    assignin('base','ShRes',ShRes);
    setflag('shading',1,flagmsg);

    fsave = get(GUI.menu.save,'Callback'); fsave();     % auto-save project
    
    if exist('cleaner','var'), cleaner(); end % if all went well, delete any partial results
end

function [ShRes,msg] = lookforprevious(Trck,SR)
% Look for existing Shading-Results (project | layout-specific)

    fmt{1} = shadingfilename('*',[],true);         % using prjpath/prjname
    fmt{2} = shadingfilename('*',Trck.name,true);  % using layout name
    [wd,~,ext] = fileparts(fmt{1});
    fmt(3:4) = strrep(fmt,ext,[ext '~']);
    if ~isequal(wd,pwd())
       fmt(5:8) = strrep(fmt,wd,'.');
    end
    candidates = cellfun(@dir,fmt,'unif',0);
    candidates = uniquecell(cat(1,candidates{:}));
    if ~isempty(candidates)
    % Refine by searching for valid TYPE keys: xGy, xIy, nP
        candidates = arrayfun(@(c) fullfile(c.folder,c.name),candidates,'unif',0);
        [~,keys,~] = cellfun(@fileparts,candidates,'unif',0);
        keys = regexpi(keys,'.*?_(\d*)([GIP]{1})(\d*)?$','tokens');
        candidates(cellfun(@isempty,keys)) = [];
    end
    
    if ~isempty(candidates)      
        setflag('shading',-3,'Checking precalculated shading results...');
        [ShRes,~,msg] = checkexistingresults(candidates,Trck,SR);
    else
        fprintf(['\tNo precalculated shading results found, (rename existing file: ',...
            shadingfilename('#P') ' or ' shadingfilename('#G/I#',Trck.name),...
            ', according to ShRes.info.type)\n']);
        msg = 'No precalculated shading results found';
        ShRes = [];
    end
end

function [ShRes,file,out] = checkexistingresults(candidates,Trck,SR)
% Evaluate 
% Minimal consistency checks for saved ShadingResults. Reduce global shading results to analysed.
    
    global SimOptions
    ShRes = [];
    
    % Check compatibility with ShadingAnalysis
    BACKUP_VARS = {'sunaz','sunel','Trackers','Terrain','ShRegions','options',...
        'simtimer','BLPoly','BshF','Nsb','DshF','DwF0','DwF','belowhorizon','t'};
    
    % Refine by searching for valid TYPE keys: xGy, xIy, nP (see ANALYSISTYPE, below)
    [~,keys,~] = cellfun(@fileparts,candidates,'unif',0);
    keys = regexpi(keys,'.*?_(\d*)([GIP]{1})(\d*)?$','tokens');
    
    % Sort candidates by type (Gridded > Irregular > Particular) and ascending resolution
    % [attempt to use high-resolution full-year results first]
    while iscell(keys{1}), keys = cat(1,keys{:}); end
    keys(cellfun(@isempty,keys)) = {'0'};
    keys = {char(keys{:,2}),...
            cellfun(@(u,d) str2double(u) + str2double(d)/10^numel(d),keys(:,1),keys(:,3))};
    keys = struct2table(struct('type',keys{1},'res',keys{2}));                            
    [keys,idx] = sortrows(keys,1:2);
    candidates = candidates(idx);
    
    if isempty(candidates)
        out = 'No valid candidates';
        return;
    end
    out = 'All candidates rejected';
    
    for j = 1:numel(candidates)
        file = candidates{j};
        
        fprintf('Checking precalculated shading results: %s\n',relativepath(file));
            
        vars = whos('-file',file);
        isbackup = isempty(setdiff(BACKUP_VARS,{vars.name}));
        if isbackup
            ShRes = load(file,'-mat');
        else
            isvalid = any(strcmp({vars.name},'ShRes'));
            if ~isvalid
                fprintf('Missing variable ShRes\n');
                ShRes = [];
                continue; 
            end
            load(file,'-mat','ShRes');
            if ~isa(ShRes,'ShadingResults')
                try
                    ShRes = ShadingResults(ShRes);
                catch ERR
                    fprintf('Failed to get SHADINGRESULTS object: %s',getReport(ERR));
                    ShRes = [];
                    continue;
                end
            end
        end

        if isbackup
            if ~isequal(ShRes.Trackers,Trck)
                fprintf('Non matching Layout');
                ShRes = [];
                continue;
            end
            if ~isequal(ShRes.ShRegions,SR)
                fprintf('Non matching ShadingRegions');
                ShRes = [];
                continue;
            end
        else
            Ntr = size(Trck.centers,2);
            Nu = numel(Trck.analysedtrackers);
            if ~any(ShRes.Nu == [Ntr,Nu])
                fprintf('ShRes.Nu (%d) does not match Trackers (%d / %d)',ShRes.Nu,Ntr,Nu);
                ShRes = [];
                continue;
            elseif Nu < Ntr && ShRes.Nu == Ntr
            % Reduce global shading results to analysed mounts
                ShRes = mountfilter(ShRes,Trck.analysedtrackers);
            end
            Np = size(Trck.analysedpoints,2);
            if ShRes.Np ~= Np
                fprintf('ShRes.Np (%d) does not match Trackers.analysedpoints (%d)',ShRes.Np,Np);
                ShRes = [];
                continue;
            end
            if ~isequaltol(ShRes.worldgeom,SR)
                fprintf('Inconsistent ShadingRegions: \n\t%s vs \n\t%s\n',...
                    ShRes.worldgeom.description,SR.description);
                ShRes = [];
                continue;
            end
            if ~isequaltol(ShRes.worldgeom,SR)
                fprintf('Inconsistent ShadingRegions: \n\t%s vs \n\t%s\n',...
                    ShRes.worldgeom.description,SR.description);
                ShRes = [];
                continue;
            end
        end

        % See ANALYSISTYPE (below) for key conventions
        switch keys{j,'type'}
        case 'G'
        % For full-year gridded results, check that resolution is compatible with SimOptions
            ds = SimOptions.shadingresolution;
            msg = sprintf('gridded shading, %0.1f° resolution',keys{j,'res'});
            if round(ds,2) ~= keys{j,'res'}
                msg = sprintf('%s (requested %0.1f°)!',msg,ds);
                warning('%s',msg);
            else
                fprintf('%s\n',msg);
            end
            if round(ds,2) < keys{j,'res'}, defanswer = 'No'; else, defanswer = 'Yes'; end
        otherwise
            msg = 'partial/point-specific shading. Might not match Meteo-Data!';
            warning(msg);
            defanswer = 'Yes';
        end    
        
        msg = sprintf('%s seems compatible: %s',relativepath(file),msg);
        if isbackup
            opt = ShRes.options;
            ShRes.options.backup = file;
            msg = [msg ', do you want to resume this calculation?']; %#ok<AGROW>
            out = ['resume from backup: ' relativepath(file)];
        else
            opt = ShRes.info.options;
            msg = [msg ', do you want to use them?']; %#ok<AGROW>
            out = ['use precalculated shading results: ' relativepath(file)];
        end

        % Compare current simulation options with those of the file
        [~,~,optmsg] = comparestruct(SimOptions,opt,'and',[],{'SimOptions','File'});
        if ~isempty(optmsg)
            warning(['The following conflicts exist between the current SimOptions and the',...
            ' results/backup in %s: \n\n%s\n'],relativepath(file),strjoin(optmsg,newline()));
        end

        switch optquestdlg(msg,'GUIshading',defanswer)
        case 'Yes'
            out = ['Attempting to ' out]; %#ok<AGROW>
            fprintf('%s\n',out);
            break;
        case 'No'
            fprintf('Rejected to %s\n',out);
            ShRes = [];
        otherwise
            error('GUIshading:questdlg','You might want to try ''Yes'' or ''No'' next time')
        end
    end
end

function resfilename = shadingfilename(key,prefix,fullpath,EXT)
% Provides consistent naming scheme for all internal subfunctions
% check compatibility with SPLITGUI.SAVESHDRES

    if nargin < 3, fullpath = false; end
    if nargin < 4, EXT = '.shdres';  end

    [prjpath,prjname] = fileparts(getSimOption('prjname'));
    if nargin < 2 || isempty(prefix), prefix = prjname; end
    resfilename = [prefix '_' key EXT];
    if fullpath, resfilename = fullfile(prjpath,resfilename); end
end

function [SunPos,type,msg] = analysistype(OPT)
% Decide on whether to do full-year(interpolation) or a partial (point-wise) analysis, based on
% the number of points for evaluation, the current settings (OPT), and possibly UI confirmation.
% Return a list of solar positions SUNPOS, standardized key TYPE to be later recognized by 
% CHECKEXISTINGRESULTS, and a feedback MSG
%
% See top GUISHADING documentation for conventions on TYPE and RESFILENAME

    global GUI
    
    ds = OPT.shadingresolution;
    N = 4*pi*sind(23.5)/(sqrt(3)*(ds*pi/180)^2); % estimated min N° vertices for a whole-year mesh

    knownlocation =  evalin('base','exist(''Location'')');
    knownsunpos = GUI.meteo.flag > 0;

    if ~knownlocation && ~knownsunpos
        msg = {'Location and Solar-Position are unknown',...
               'Import Meteo-Data, or provide a *.tif DEM file as Terrain'};
        setflag('shading',-1,msg);
        error(strjoin(msg,newline()));
    end

    if knownsunpos
        SunPos = evalin('base','SunPos');
        minSunEl = OPT.minSunEl;
        minGHI = OPT.minGHI;
        GHI = evalin('base','MD.GHI');
        notdark = SunPos.El >= minSunEl;
        trouble = (GHI > minGHI) & ~notdark;
        if any(trouble)
            warning(['%d time-steps below minSunEl (%0.1f°) with GHI > threshold (%0.1f W). ',...
                     'GHI values range from %0.1f to %0.1f W, will yield GTI = 0.'],...
                     nnz(trouble),minSunEl,minGHI,min(GHI(trouble)),max(GHI(trouble)));
        end
        notdark = notdark & (GHI > minGHI);
        gridded = nnz(notdark) > 0.8*N;   
        if ~gridded && knownlocation
            switch optquestdlg(sprintf(['Full-year (interpolation) shading analysis is expected ',...
                'to require %0.1fk points, whereas meteo-data has only %d daylight points. ',...
                'Do you want to perform a (faster) partial shading-analysis? the results in ',...
                'this case would be meteo-data-specific'],N/1000,nnz(notdark)),...
                'GUIshading','Full','Partial','Partial')
            case 'Full', gridded = true;         % force use of grid
            case 'Partial', knownlocation = false;
            otherwise, error('Stopped by user');
            end
        end
    else
        gridded = true; % force use of grid
    end

    if gridded      
        if knownsunpos, Nt = numel(SunPos.Az); end
        gridspecs = sprintf('%0.1f° (%0.1f min)',ds,ds*4);        
        if knownlocation
        % Replace any existing solar positions with a uniform full-year grid
            msg =  ['Using full-year grid at ',gridspecs];
            Loc = evalin('base','Location');
            [az,el] = solarposition.sunposgrid(Loc.latitude,ds,'minSunEl',OPT.minSunEl);
            SunPos = struct('Az',az,'El',el);
            
            type = strrep(num2str(ds,'%0.3g'),'.','G'); % 1G25 = 1.25° Global Grid
        else
        % Generate an interpolation mesh from existing positions
            msg = [gridspecs, ' mesh from scattered solar-positions (unkown location)'];
            notdark = SunPos.El >= OPT.minSunEl;
            SunPos.Az = SunPos.Az(notdark); 
            SunPos.El = SunPos.El(notdark);
            V = sunposmesh(SunPos.Az,SunPos.El,'meshsize',(ds*pi/180),'offset',(ds*pi/180)/2);
            az = atan2d(V(:,2),V(:,1));
            el = atan2d(V(:,3),hypot(V(:,1),V(:,2)));
            SunPos = struct('Az',az,'El',el);
            
            type = strrep(num2str(ds,'%0.3g'),'.','I'); % 3I5 = 3.5° Irregular/Incomplete mesh
        end

        if knownsunpos
            msg = sprintf('%s: %d mesh vertices, from %d original points.',...
                msg,numel(SunPos.Az),Nt);
        else
            msg = sprintf('%s: %d mesh vertices.',msg,numel(SunPos.Az));
        end
    else
    % ... or, keep the original solar positions
        type = sprintf('%dP',numel(SunPos.Az)); % 4000P = 4000 Particular Points
        msg = [type,' Point-specific shading analysis.'];
    end
    fprintf('\n%s\n',msg);
end