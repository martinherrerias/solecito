function varargout = splitscript(prjname,varargin)
% SPLITSCRIPT(PRJNAME,..) - automated, consecutive run of simulation steps, to be run  
% instead of SPLITGUI for non-user-interactive (e.g. batch-server) applications. The script is 
% intended to be run in -nodesktop, -noFigureWindows mode, typing from the System console:
%
%   matlab -nosplash -nodesktop -noFigureWindows -r splitscript(...)
%
%   SPLITSCRIPT(PRJNAME) will open an existing *.ssp project file PRJNAME (or create a new one
%   with that name if it doesn't exist), and sequentially try to excecute all simulation steps 
%   whose flag is set to zero (Ready-to-Run).
%   For a new project, this means starting from 'setup' and trying to go all the way to 'solver'
%   For a half-completed project, excecution will resume from the first incomplete task.
%
%   NOTE: in this mode, SPLITSCRIPT will not re-run steps that have already been completed.
%
% SPLITSCRIPT() will run any and all *.ssp project files in PWD sequentially, or, if no existing
%   create an auto-named new project in the working directory, then try to run a simulation from 
%   start to finish.
%
% SPLITSCRIPT(PRJNAME,STEPS) where STEPS is a cell-array of labels, will try to sequentially run
%   all steps listed in STEPS, REGARDLESS of their status flag (!). Step labels are:
%
%            'setup' : Options
%         'phystrck' : Mount - Physical
%          'physmod' : Module - Physical
%           'layout' : Layout
%          'terrain' : Terrain
%           'arrdef' : Array Definition
%           'models' : Electrical Models
%            'meteo' : Meteo-data
%          'shading' : Shading-Analysis
%         'irrtrans' : POA-Irradiance (§)
%           'solver' : Circuit Solution
%
%   NOTE: steps will try to be executed in the order in which they are provided, and execution 
%   will stop upon occurence of the first error (see 'onerror' below).
%
% SPLITSCRIPT(PRJNAME,FROMSTEP,TOSTEP) works similarly, running all(§) steps between the labels
%   FROMSTEP and TOSTEP, regardless of their status flag.
%
%   (§) An exception is the memory intensive POA-irradiance step 'irrtrans', which will by default
%   be skipped to be performed internally by 'solver' (that is, whenever 'solver' is in the queue).
%   See documentation for GUISOLVER and GUIIRRTRANS for details.
%
% SPLITSCRIPT(..,'onerror',ACTION) - 'debug', 'continue', or 'stop' on error. The default option 
%   is 'continue', i.e. attempt all remaining steps, even if some fail. If dbstop-on-error is set
%   the default changes to 'debug', to reflect this behavior. 'stop' halts execution after the 
%   first error.
%
% SPLITSCRIPT(..,'simoptions',{'name',val,...}) - Use SETSIMOPTION('name',val) for every, name-
%   value pair provided, prior to executing any other simulation step.
% 
% See also: SPLITGUI, SETFLAG

    global GUI
    
    % Parse project name
    if nargin < 1 || isempty(prjname)
        prjname = pickfile('*.ssp',Inf);
        if isempty(prjname), prjname = ['newproj_' datestr(now(),'yymmdd_HHMMSS') '.ssp']; end
    else
        prj = pickfile(prjname,Inf);
        if ~isempty(prj)
            prjname = prj;
        end
    end 
    if iscell(prjname) && isscalar(prjname), prjname = prjname{1}; end
    if iscell(prjname)
        varargout{1} = cellfun(@(p) splitscript(p,varargin{:}),prjname);
        return;
    end
    assert(ischar(prjname) && strcmp(prjname(end-3:end),'.ssp'),...
        'splitscript:prjname','Project name must be a string, with extension *.ssp');
    
    if ~isempty(dbstatus), opt.onerror = 'debug'; else, opt.onerror = 'continue'; end
    opt.simoptions = {};
    [opt,varargin] = getpairedoptions(varargin,opt);

    switch lower(opt.onerror)
        case 'stop', stoponerror = true;
        case 'continue', stoponerror = false;
        case 'debug'
            dbstop if error;
            stoponerror = false;
        otherwise
            error('Unrecognized ''onerror'' option.');
    end
    
    % Backup base workspace and reset at the end (splitGUI will clear it)
    basedir = pwd(); 
    if ~evalin('base','isempty(whos())')
        evalin('base','save(''.splitscript_backup.mat'');');
        lastwill = onCleanup(@() reloadworkspace(basedir));
    else
        lastwill = onCleanup(@() cd(basedir));
    end
    
    splitGUI(prjname);
    fprintf('Running from SPLITSCRIPT(''%s'')\n',prjname);
    
    % parse step-list
    allsteps = setdiff(fieldnames(GUI),{'menu'},'stable');
    switch numel(varargin)
        case 0
            steps = allsteps;
            unknown = {};
            forcerun = false;
            explicit = false;
        case 1
            assert(iscellstr(varargin{1}),'splitscript:cellarg',...
                'In splitscript(PRJNAME,STEPS), STEPS must be a cell-array of strings');
            [steps,unknown] = match(varargin{1},allsteps); 
            forcerun = true;
            explicit = true;
        case 2
            if isempty(varargin{1}), varargin(1) = allsteps(1); end
            if isempty(varargin{2}), varargin(2) = allsteps(end); end
            assert(ischar(varargin{1}) && ischar(varargin{2}),'splitscript:charargs',...
                'In splitscript(PRJNAME,FROM,TO) all arguments must be strings');
            [steps,unknown,n] = match(varargin(1:2),allsteps);
            if any(n == 0)
                % leave as it is, error will be reported below (*)
            elseif n(2) < n(1)
                error('splitscript:tofrom','''%s'' is after ''%s'', check argument order',steps{:});
            else
                steps = allsteps(n(1):n(2));
            end
            forcerun = true;
            explicit = false;
        otherwise
            error('Unrecognized arguments')
    end
    assert(isempty(unknown),'splitscript:steps',shortliststr(unknown,'Unknown step',5,'quotes',''''));
    
    % Use in-solver transposition, unless requested
    if ~explicit && ismember('solver',steps)
        steps = setdiff(steps,{'irrtrans'},'stable');
    end
    
    setSimOption('version','batch');
    setSimOption('plotting',false);
    
    % Override SimOptions if provided as 'simoptions',{'opt1',val,'opt2',val2,...}
    if ~isempty(opt.simoptions)
        try opt.simoptions = reshape(opt.simoptions,2,[]); catch, opt.simoptions = NaN; end
        assert(iscellstr(opt.simoptions(1,:)),'Non-matching SimOption pairs'); %#ok<ISCLSTR>
        cellfun(@setSimOption,opt.simoptions(1,:),opt.simoptions(2,:));
    end
    
    % Run simulation steps
    flags = nan(numel(steps),1);
    for j = 1:numel(steps)
    % for each listed step x, run GUIx(), e.g. GUIsetup(), GUIphysmod, etc.
        if GUI.(steps{j}).flag > 0 && ~forcerun, continue; end
        
        ERR = MException.empty;
        if strcmpi(opt.onerror,'debug')
        % see what happens (dbstop if error)
            eval(['GUI' steps{j} '();']);
        else
        % exit/continue on error
            try eval(['GUI' steps{j} '();']);
            catch ERR
                if ~stoponerror, disp(getReport(ERR)); end
            end
        end
        flags(j) = GUI.(steps{j}).flag;
        if ~isempty(ERR), flags(j) = -1; end          % should never come to this, but hey...
        if flags(j) == -1 && stoponerror, break; end  % stop at the first error
        
        if isequal(steps{j},'setup') && ~isempty(opt.simoptions)
        % Override options (again) after setup
            cellfun(@setSimOption,opt.simoptions(1,:),opt.simoptions(2,:));
        end
    end
        
    % Save project, whether it crashed or finished smoothly
    fsave = get(GUI.menu.save,'Callback');
    fsave();
    % close('all');
    close(GUIfigure('GUImain'));
    
    if nargout > 0, varargout{1} = any(flags == -1)*1; end  % exit flag = 0: ok, 1: error

    if flags(j) == -1 && stoponerror
        rethrow(ERR);
    end
    
    % clear lastwill; % restore original workspace
end

function [in_B,not_in_B,idx_B] = match(a,b)
    [ia,idx_B] = ismember(lower(a),lower(b));
    in_B = b(idx_B(ia));
    not_in_B = a(~ia);
end

function reloadworkspace(basedir)
    cd(basedir);
    evalin('base','clear variables'); 
    evalin('base','load(''.splitscript_backup.mat'');');
    delete('.splitscript_backup.mat');
end
