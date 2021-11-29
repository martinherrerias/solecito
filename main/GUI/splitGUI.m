function splitGUI(prjname)
% SPLITGUI() - Main simulation script, minimal GUI to divide simulation steps into independent but
% related tasks. Each step is represented by single button (that links to a dedicated function)
% and a text box for short output. The structure of SPLITGUI is designed so that the number of 
% simulation steps, their order and mutual relations are easily customized using internal functions
% MODSTEP and ADDREQ. Although the function doesn't directly return a result, it does create the
% global variables GUI and SimOptions. For each step 'stepx', GUI has the following structure:
%
%   GUI.stepx.btn - handle for button object for stepx
%   GUI.stepx.txt - handle for text object for stepx
%   GUI.stepx.flag - status flag for stepx (see SETFLAG):     -3: Running (disabled)
%                                                             -2: Disabled
%                                                             -1: Error in last run (enabled)
%                                                              0: Ready to run
%                                                              1: Completed without errors
%                                                              2: Completed with warnings
%
%   GUI.stepx.req - vector of indices, pointing at the required steps that need to be completed
%       (flag > 0), in order for stepx to be enabled.
%
%   global SimOptions is a structure with field {version} which after GUIsetup will hold
%   all simulation options (precision, simplifications, models, etc., see GUISETUP).
%
% Functions associated with each step are expected to load any variables they need from base 
% (or external files), perform operations inside a try-catch block (so that the crash of one stage
% doesn't require a new complete simulation), and write results to the base-workspace, or external
% files.
%
% SPLITGUI(PRJNAME) automatically loads/creates a project PRJNAME (path and filename)
%
% See also: SETFLAG, UPDATEFLAGS, GUISETUP, GUIPHYSTRCK, GUIPHYSMOD, GUILAYOUT, GUITERRAIN,
% GUIARRDEF, GUIREDUCE, GUIMETEO, GUISHADING, GUIIRRTRANS, GUIMODELS, GUISOLVER

    feature('DefaultCharacterSet','UTF-8');
    
    evalin('base','clear variables');
    clearvars -global -except GUI
    clear PWLparallel PWLseries mkivpp parsetolerance % clear persistent
    clear stopwatch radiance_igawa onediodeMPP onediodemodel onediodemodel2
    
    warning on all
    % warning off backtrace
    % warning off verbose
    
    global SimOptions;

    SimOptions = DefaultOptions();      % provisional (*)
    SimOptions.version = 'UI';
    %clc
    
    % Record everything on log file
    setdiary(); % start on provisional diary file
    fprintf('\nSplitGUI: Log started on the %s\n', datestr(now,'dd.mm.yyyy-HH:MM:SS'));
 
    figh = GUIfigure('GUImain','SplitGUI','main');
    drawform(figh); % Draw GUI, details below
    
    SimOptions = struct('version',SimOptions.version); %(*)
    
    if nargin > 0 && ~isempty(prjname)
        if isempty(dir(prjname))
            newprj([],[],prjname);
        else
            prjname = pickfile(prjname);
            f = fopen(prjname,'r'); 
            assert(f > 0,'splitGUI','Could not open project %s',prjname);
            fclose(f);
            openprj([],[],prjname);
        end
    end
end

function drawform(figh)
% Draw a form with two columns, each row composed of a button btn(i) and a textbox stat(i)
% Associate functions to buttons, set starting captions, and define relationships (req) among steps

    global GUI
    GUI = struct();
    evalin('base','global GUI');
    
    figure(figh); clf(figh);
    set(figh,'CloseRequestFcn',@closeform);
    
    set(figh,'HandleVisibility','on'); % provisionally let elements be added
    
    % Create menu
    set(figh,'MenuBar','none');
    GUI.menu.proj = uimenu('Label','Project');
    GUI.menu.new = uimenu(GUI.menu.proj,'Label','New','Accelerator','N','Callback',@newprj);
    GUI.menu.open = uimenu(GUI.menu.proj,'Label','Open','Accelerator','O','Callback',@openprj);
    GUI.menu.save = uimenu(GUI.menu.proj,'Label','Save','Separator','on','Accelerator','S','Callback',@saveprj,'Enable','off');
    GUI.menu.saveas = uimenu(GUI.menu.proj,'Label','Save As...','Callback',@(~,~) saveprjas(),'Enable','off');
    GUI.menu.minprj = uimenu(GUI.menu.proj,'Label','Save minimal...','Callback',@minprj,'Enable','off');
    GUI.menu.reduce = uimenu(GUI.menu.proj,'Label','Save reduced...','Callback',@(~,~) GUIreduce(),'Enable','off');
    GUI.menu.quit = uimenu(GUI.menu.proj,'Label','Quit','Separator','on','Accelerator','Q','Callback',@quitgui);
    
    GUI.menu.tools = uimenu('Label','Tools');
    GUI.menu.simopt = uimenu(GUI.menu.tools,'Label','Reload Options','Callback',@GUIsetup,'Enable','off');
    GUI.menu.export = uimenu(GUI.menu.tools,'Label','Export...','Separator','on','Enable','on');
    GUI.menu.export_shdres = uimenu(GUI.menu.export,'Label','Shading Results...','Callback',@saveshdres,'Enable','off');
    GUI.menu.export_solres = uimenu(GUI.menu.export,'Label','Solver Results...','Callback',@savesolres,'Enable','off');
    GUI.menu.export_balance = uimenu(GUI.menu.export,'Label','Energy Balance...','Callback',@savebalance,'Enable','off');
    GUI.menu.print_summary = uimenu(GUI.menu.tools,'Label','Result Summary','Separator','on','Callback',@GUIprintsummary,'Enable','off');
    
    % GUI.menu.parallel = uimenu('Label','Parallel');
    % GUI.menu.divide = uimenu(GUI.menu.parallel,'Label','Divide Project...','Callback',@GUIdivide,'Enable','off');
    % GUI.menu.runparallel = uimenu(GUI.menu.parallel,'Label','Solve parallel...','Callback',@GUIrunparallel,'Enable','off');
    % GUI.menu.merge = uimenu(GUI.menu.parallel,'Label','Merge Results...','Callback',@GUImerge,'Enable','off');
    
    GUI.menu.help = uimenu('Label','Help');
    GUI.menu.doc = uimenu(GUI.menu.help,'Label','Documentation...','Accelerator','H','Callback','','Enable','off');
    GUI.menu.about = uimenu(GUI.menu.help,'Label','About...','Separator','on','Callback','','Enable','off');

    btndims = [1,1];                   % button width, height
    textdims = [3,btndims(2)];         % texbox width, [height]
    marg = [0.2,0.2,0.4,0.4];          % [l r t b] margins
    gap = [0.2,0.1];                   % gap between rows and columns
    
    idlist = cell(0); 
    labels = cell(0);
    btnfuncts = cell(0);
    flags = [];
    req = cell(0);
    
    % Define simulation steps
    modstep('setup','Options',{'Create a new project, or load an existing *.ssp file',...
                               'Use an AddSimOptions.m file for custom options.'});
    modstep('phystrck','Mount - Physical','*.tpoly files required.');
    modstep('physmod','Module - Physical','*.mpoly or *.modint file required.');
    modstep('layout','Layout',' *.mounts file required.');
    modstep('terrain','Terrain','Uses optional *.hor/.*tif file(s)');
    modstep('arrdef','Array Definition','*.def file required.');
    modstep('models','Electrical Models','*.samlib files required.');
    modstep('meteo','Meteo-data','Requires *.meteo file');
    modstep('shading','Shading-Analysis','Requires complete layout');
    modstep('irrtrans','POA-Irradiance','Requieres Shading Analysis & Meteo Data');
    modstep('solver','Circuit Solution','All steps above required.');
   
    % Define requirements for each step
    addreq('setup',{});
    addreq('phystrck',{'setup'});
    addreq('physmod',{'phystrck'});
    addreq('layout',{'phystrck'});
    addreq('terrain',{'layout'});
    addreq('arrdef',{'layout'});
    addreq('models',{'setup'});
    addreq('meteo',{'setup'});
    addreq('shading',{'physmod','terrain'});
    addreq('irrtrans',{'shading','meteo'});
    addreq('solver',{'shading','meteo','models'});

    Nb = numel(idlist);     % number of rows
    
    % Get element positions
    h = (marg(3)+marg(4)+Nb*btndims(2)+(Nb-1)*gap(2));
    w = (marg(1)+marg(2)+btndims(1)+textdims(1)+gap(1));
    xx = marg(1)/w + [0,btndims(1)+gap(1)]/w;
    yy = 1-(marg(3) + (1:Nb)*btndims(2)+(0:Nb-1)*gap(2))/h;

    % Fill-up buttons, labels, and GUI structure
    for j = 1:Nb
        pos = [xx(1) yy(j) btndims./[w,h]];
        %if flags(j)>= 0, enab = 'on'; else, enab = 'off'; end
            
        btn = uicontrol('Parent',figh,'Style','pushbutton','String',labels{j,1},...
            'Units','normalized','Position',pos,'Callback',btnfuncts{j},'Enable','off');

        pos = [xx(2) yy(j) textdims./[w,h]];
        stat = uicontrol('Parent',figh,'Style','text','String',labels{j,2},...
            'Units','normalized','Position',pos,'Visible','on','HorizontalAlignment','left');

        GUI.(idlist{j}).btn = btn;
        GUI.(idlist{j}).txt = stat;
        GUI.(idlist{j}).lbl0 = labels{j,2};
        GUI.(idlist{j}).req = req{j};
        setflag(idlist{j},flags(j),labels{j,2});
    end
    
    set(figh,'HandleVisibility','off'); % block as target for graphics output
    
    updateflags() % enable any steps that have no requirements, list requirements for the rest
            
    function modstep(id,btnlabel,txtlabel)
    % Add a simulation step (a button and textbox in the GUI)
    % Button Callback function names are @GUI<id> (as assumed by RESTOREHANDLES)
    
        btnfunc = str2func(['@GUI' id]);
        idx = find(strcmp(id,idlist));
        if isempty(idx), idx = numel(idlist)+1; end

        idlist{idx} = id;
        labels{idx,1} = btnlabel;
        labels{idx,2} = txtlabel;
        btnfuncts{idx} = btnfunc;
        flags(idx)= -2; % start with all steps disabled
        req{idx} = [];
    end

    function addreq(id,reqlist)
    % Add a list of requirements (other simulation steps) to step ID
    
        idx = strcmp(id,idlist); 
        rq = cellfun(@(s) find(strcmp(s,idlist)),reqlist);
        req{idx} = rq(:)';
    end
end

function closeform(~,~)
    
    fclose('all');
    delete(findall(0,'Tag','TMWWaitbar'));
    
    % Close all associated figure windows
    warning off 'getSimOption:getdefault'
    % fh = getSimOption('fighandles');
    fh = GUIfigure('-all');
    for f = fieldnames(fh)'
        if strcmp(f{1},'GUImain'), continue; end
        if ishandle(fh.(f{1})), close(fh.(f{1})); end 
    end
    
    fprintf('\nSplitGUI: Log ended on the %s\n', datestr(now,'dd.mm.yyyy-HH:MM:SS'));
    diary OFF
    
    closereq();
end

function quitgui(~,~)
    close(GUIfigure('GUImain'))
end

function saveprj(~,~)
    global GUI
    
    allsteps = setdiff(fieldnames(GUI),{'menu'},'stable');
    
    if isempty(getSimOption('prjname')),saveprjas(); return; end
    
    % Replace txt-handles by txt-strings in saved workspace copy
    for s = allsteps', GUI.(s{1}).txt = get(GUI.(s{1}).txt,'String'); end
    lastwill = onCleanup(@restorehandles); % restore txt-handles when done
    
        evalin('base','save(getSimOption(''prjname''));'); % save workspace
    
    fprintf('\nProject saved to: %s\n\n',getSimOption('prjname'));
end

function saveprjas(filename)

    global GUI
    allsteps = setdiff(fieldnames(GUI),{'menu'},'stable');
    previousname = getSimOption('prjname');
    
    if nargin < 1 || isempty(filename)
        if isempty(previousname), filename = '*.ssp';
        else
            [path,filename,ext] = fileparts(previousname);
            filename = fullfile(path,[filename '_copy' ext]);
        end
    end
    [filename,path] = uiputfile(filename);
    if isequal(filename,0), return; end

    filename = strrep(fullfile(path,filename),'\','/'); % avoid scape-char issues
    setSimOption('prjname',filename);
    set(GUIfigure('GUImain'),'Name',getSimOption('prjname'));
       
    % Replace txt-handles by txt-strings in saved workspace copy
    for s = allsteps', GUI.(s{1}).txt = get(GUI.(s{1}).txt,'String'); end
    lastwill = onCleanup(@restorehandles); % restore txt-handles when done
    
    evalin('base','save(getSimOption(''prjname''));'); % save workspace
    
    fprintf('\nProject saved to: %s\n\n',getSimOption('prjname'));
    
    % Switch to new diary
    setdiary(filename);
    fprintf('\nSplitGUI: Log started on the %s\n', datestr(now,'dd.mm.yyyy-HH:MM:SS'));
    fprintf('\nProject saved-as from: %s\n\n',previousname);
end

function minprj(~,~)
% Save a minimal copy of the project, ready for shading...solver

    global GUI
    
    % FUTURE: the list should be derived from a GUI.varin/varout field similar to GUI.req
    USEDVARS = {'SimOptions','GUI','OptionFlags','Location',...
                'Time','SunPos','MeteoData','HorizonProfile',...
                'Trackers','ArrayDef','ArrIdx',...
                'ModIVint','Diode','Inverter','celltemp'};
    
    %saveproj();  % auto-save project
    
    % Clear everything but USEDVARS
    v = evalin('base','whos');
    v = setdiff({v.name},USEDVARS);
    evalin('base',['clear ' strjoin(v,' ')]);
    
    for simstep = {'solver','irrtrans','shading'}
        setflag(simstep{:},0,GUI.(simstep{:}).lbl0);
    end

    [filename,path] = uiputfile('*.ssp');
    if isequal(filename,0), return; end
    
    previousname = getSimOption('prjname');
    
    filename = strrep(fullfile(path,filename),'\','/'); % avoid scape-char issues
    setSimOption('prjname',filename);
    set(GUIfigure('GUImain'),'Name',getSimOption('prjname'));
       
    % Replace txt-handles by txt-strings in saved workspace copy
    allsteps = setdiff(fieldnames(GUI),{'menu'},'stable');
    for s = allsteps', GUI.(s{1}).txt = get(GUI.(s{1}).txt,'String'); end
    lastwill = onCleanup(@restorehandles); % restore txt-handles when done
    
    evalin('base','save(getSimOption(''prjname''));'); % save workspace
        
    fprintf('\nProject saved to: %s\n\n',getSimOption('prjname'));
    
    % Switch to new diary
    setdiary(filename);
    fprintf('\nSplitGUI: Log started on the %s\n', datestr(now,'dd.mm.yyyy-HH:MM:SS'));
    fprintf('\nMinimal copy of: %s\n\n',previousname);
end

function openprj(~,~,filename)
    global GUI
    global SimOptions
    global OptionFlags
    
    PERSISTENT = {'version'}; % session-dependent variables
    
    % Save session-dependent variables in backup options structure...
    opt = cat(1,PERSISTENT,cellfun(@getSimOption,PERSISTENT,'unif',0));
    opt = struct(opt{:});
    
    GUI_copy = GUI; % save a copy, to make persistent
    
    allsteps = setdiff(fieldnames(GUI),{'menu'},'stable');
    
    if nargin < 3 || isempty(filename)
        filename = pickfile('*.ssp','ui',true,'fullpath',true);
    else
        filename = pickfile(filename,'fullpath',true);
    end
    filename = strrep(filename,'\','/'); % avoid scape-char issues
    cd(fileparts(filename));
    
    % ** RESTORE WORKSPACE **
    evalin('base','clearvars'); 
    clear PWLparallel PWLseries mkivpp parsetolerance % clear persistent
    
    % Load/complete SimOptions & OptionFlags first (version compatibility): 
    load(filename,'-mat','SimOptions','OptionFlags');
    evalin('base','global GUI SimOptions OptionFlags;');
    completeoptions(); % update/complete options structure to current

    % Load (almost) everything else
    allvars = whos('-file',filename); 
    allvars = {allvars.name};
    vars = setdiff(allvars,{'SimOptions','OptionFlags','ans','ShRes'});
    evalin('base',sprintf('load(''%s'',''-mat'',''%s'');',filename,strjoin(vars,''',''')));
    
    % Change from PVL_MAKETIMESTRUCT to DATETIME (2020.10.20)
    if any(strcmp(allvars,'Time')) && any(strcmp(allvars,'MeteoData')) && ...
            ~evalin('base','isdatetime(Time)')
        T = evalin('base','Time');
        MD = evalin('base','MeteoData');
        if isfield(T,'dt')
            dt = T.dt;
            T = parsetime(T,'step',dt);
        else
            [T,dt] = parsetime(T);
        end
        MD.timestep = dt;
        assignin('base','Time',T);
        assignin('base','MeteoData',MD);
    end
    
    % Load shading-results last (allows update of deprecated versions)
    if any(strcmp(allvars,'ShRes'))
        evalin('base',sprintf('load(''%s'',''-mat'',''ShRes'');',filename));
    end
    % **

    % Copy session-dependent options into SimOptions
    cellfun(@(s) setSimOption(s,opt.(s)),PERSISTENT);

    % Copy loaded txt-labels and flags from loaded project. Dump the rest (use session GUI)
    cellfun(@(s) set(GUI_copy.(s).txt,'String',GUI.(s).txt),allsteps);
    for s = allsteps(:)', GUI_copy.(s{1}).flag = GUI.(s{1}).flag; end
    GUI = GUI_copy;
    
    setdiary(filename); % update diary location
    
    updateflags();
    setSimOption('prjname',filename); % in case either has changed
    fprintf('\n%s - Project loaded from: %s\n\n',datestr(now,'dd.mm.yyyy-HH:MM:SS'),filename);
    set(GUIfigure('GUImain'),'Name',relativepath(filename));    
end

function newprj(~,~,filename)
    global GUI
    global SimOptions
    %
    allsteps = setdiff(fieldnames(GUI),{'menu'},'stable');
    
    nameisvalid = @(s) isempty(regexp(s,'[*?"<>|]','once'));
    if nargin < 3 || isempty(filename) || ~nameisvalid(filename)
        [filename,path] = uiputfile('*.ssp');
        if isequal(filename,0), return; end
    else
        path = './';
    end
    filename = fullfile(path,filename);
    filename = strrep(filename,'\','/'); % avoid scape-char issues
    
    % Restart GUI unless just starting... 
    if any(cellfun(@(s) GUI.(s).flag,allsteps)~= -2)
        splitGUI();
        GUI = evalin('base','GUI');
    end
    
    SimOptions.prjname = filename;
    setdiary(filename); % update diary location
    GUIsetup(); % Load options
    
    saveprj();    
end

function saveshdres(~,~)
% Export ShRes as LAYOUT_TYPE.shdres, where LAYOUT or PRJNAME_nP.shdres, depending on
% whether ShRes.info.layout and ..type are available.

    PREFIX = '';       
    EXT = '.shdres';  % check compatibility with GUISHADING
    
    [prjpath,prjname] = fileparts(getSimOption('prjname'));
    
    ShRes = evalin('base','ShRes');
    
    if all(isfield(ShRes.info,{'layout','type'}))
        filename = fullfile(prjpath,[PREFIX ShRes.info.layout '_' ShRes.info.type EXT]);
    else
        type = sprintf('%dP',ShRes.Nt);
        filename = fullfile(prjpath,[PREFIX prjname '_' type EXT]);
    end

    save(filename,'ShRes');
    fprintf('\nShading-Results saved to: %s\n\n',filename);
end

function savesolres(~,~)
    
    [prjfolder,prjname] = fileparts(getSimOption('prjname'));
    filename = fullfile(prjfolder,['solres_' prjname '.mat']);
    evalin('base',['save(''' filename ''',''SolRes'');']);
    fprintf('\nSolver-Results saved to: %s\n\n',filename);
end

function savebalance(~,~)
    
    [prjfolder,prjname] = fileparts(getSimOption('prjname'));
    filename = fullfile(prjfolder,['EnergyBalance_' prjname '.xls']);
    EB = evalin('base','EB');
    Nt = size(EB.GHI,1);
    struct2csv(EB,Nt,filename);
    fprintf('\nEnergy-Balance saved to: %s\n\n',filename);
end

% function GUIdivide(~,~)
% % PROVISIONAL: shouldn't be necessary when parfor loops are added
% 
%     global GUI
%     
%     p = gcp();
%     N = p.NumWorkers;
%     fprintf('\nDividing Project into Workers...\n');
% %     if runningfromUI()
% %         N = inputdlg('Input number of pieces for Project-Division','Divide...',1,{'4'});
% %         N = str2double(N);
% %         if isempty(N)||isnan(N), fprintf('\bor maybe not... -_-\n'); return; end
% %     else
% %         N = 4;
% %     end
%     
%     saveproj();
%     divideproject(getSimOption('prjname'),N);
%     fprintf('\t...Project split into %d pieces.\n',N);
%     set(GUI.menu.runparallel,'Enable','on');
% end
% 
% function jobs = GUIrunparallel(~,~,varargin)
%     global GUI
%     DIR_NAME = 'sub';
%     WBTAG = '<WAITBAR>';
%     
%     prjname = getSimOption('prjname');
%     [prjpath,prjname] = fileparts(prjname);
%     
%     % Find part-files...
%     subs = dir(fullfile(prjpath,DIR_NAME,[prjname '_*.ssp']));
%     subs([subs.isdir]) = []; % keep files only
%     subs = cellfun(@(s) fullfile(prjpath,DIR_NAME,s),{subs.name},'unif',0);
% 
%     gcp();
%     jobs = cellfun(@(s) parfeval(@splitscript,0,s,varargin{:}),subs);
%     % wait(jobs);
%     % opt = {'CaptureDiary',false};
%     % jobs = cellfun(@(s) batch(@splitscript,0,[{s},varargin,opt]),subs); 
% 
%     assignin('base','jobs',jobs);
%     finished = false;
%         
%     % logID = zeros(numel(jobs),1);
%     mergedlog = cell(0,1);
%     lastlogline = zeros(numel(jobs),1);
%     diaries = cell(numel(jobs),1);
%     newlines = cell(numel(jobs),1);
%     wblines = cell(numel(jobs),1);
%     while ~all(finished)
%         pause(4);
%         wblines(:) = {cell.empty};
%         finished = arrayfun(@(j) strcmp(j.State,'finished'),jobs);
%         diaries(finished) = {cell.empty};
%         [diaries{~finished}] = deal(jobs(~finished).Diary);
%         pause(1);
%         
%         for j = find(~cellfun(@isempty,diaries))'       
%             % if logID(j) <= 0
%             %     try logID(j) = fopen([subs{j}(1:end-4) '.log'],'r');
%             %     catch, continue;
%             %     end
%             % end
%             % data =  textscan(logID(j),'%[^\n]');
%             % newlines{j} = data{1};
% 
%             newlines{j} = strsplit(diaries{j},newline())';
%             newlines{j}(1:lastlogline(j)) = [];
%             if isempty(newlines{j}), continue; end
%             lastlogline(j) = lastlogline(j) + numel(newlines{j});
% 
%             % Keep #waitbar lines for later...
%             waitbars = contains(newlines{j},WBTAG);
%             wblines{j} = newlines{j}(waitbars);
%             newlines{j}(waitbars) = [];
%             if isempty(newlines{j}), continue; end
% 
%             % Remove lines which have already been printed...
%             % ... (or will be printed for other jobs)
%             printed = cellfun(@(s) any(strcmp(s,mergedlog)),newlines{j});
%             newlines{j}(printed) = [];
%             if isempty(newlines{j}), continue; end
%             mergedlog = cat(1,mergedlog,newlines{j});
%             
%             cellfun(@(s) fprintf('Job %02d: %s\n',j,s),newlines{j}); 
%         end
%         
%         if all(cellfun(@isempty,newlines)) && ~all(cellfun(@isempty,wblines))
%         % Print the WAITBAR line with max progress
%             allwblines = cat(1,wblines{:});
%             prog = regexp(allwblines,'<WAITBAR>\(([\d.]+)%\).*','tokens');
%             prog = cat(1,prog{:});
%             [~,idx] = max(str2double(cat(1,prog{:})));
%             fprintf('%s\n',allwblines{idx});
%         end
%     end
% 
%     set(GUI.menu.merge,'Enable','on');
% end
% 
% function GUImerge(~,~)
%     prjname = getSimOption('prjname');
%     mergeproject(prjname);
%     [prjpath,prjname] = fileparts(prjname);
%     openproj([],[],fullfile(prjpath,[prjname '_merge.ssp']));
% end

function GUIprintsummary(~,~)
    printsummary({'irrtrans','solver'});
end

function setdiary(prjname)
% SETDIARYFILE() - starts DIARY on a (provisional) ./~solarsim.log file.
% SETDIARYFILE(prjname) - switch to DIARY('./prjname.log'), if previous file was ./~solarsim.log
%       temporary file, copy its contents to new file. 

    DEFNAME = 'solarsim.log~';
    global logfile;
    if isempty(logfile),logfile = ''; end
 
    if nargin == 0
        newfile = fullfile(pwd(),DEFNAME);
        if exist(newfile,'file'), backupdelete(newfile); end % clear forgotten logs
    else
        [path,newfile,~] = fileparts(prjname);
        newfile = fullfile(path,[newfile '.log']);
    end
    
    diary off
    
    if ~isempty(regexpi(logfile,[DEFNAME '$'])) && exist(logfile,'file')
        % Copy contents of provisional log-file
        fin=fopen(logfile,'r');
        temp = fread(fin,'uint8');
        fclose(fin);

        % Append contents to existing file, or create new one
        fout = fopen(newfile,'a+'); 
        fwrite(fout,temp,'uint8');
        fclose(fout);
        delete(logfile);
        
        % Remove <WAITBAR>(..)... lines
        cleanlog(newfile);
    end
    
    diary(newfile);
    logfile = newfile;
end

% function recoverhandles(GUI)
%     allsteps = setdiff(fieldnames(GUI),{'menu'},'stable');
%     n = numel(allsteps);
%     c = flipud(get(GUIfigure('GUImain'),'children'));  % children are returned first-in-last-out
%     c = c(4:end);                     % remove menu entries
% 
%     for j = 1:n
%         GUI.(allsteps{j}).btn = c(2*(j-1)+1);
%         GUI.(allsteps{j}).txt = c(2*j);
%     end
% end
