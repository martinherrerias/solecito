function updateflags(simstep)
% UPDATEFLAGS(SIMSTEP) - Evaluates the current GUI flags (GUI.(s).flag for all steps s), and
% determines if simulation step SIMSTEP is ready to be enabled, based on GUI.(SIMSTEP).req
% If SIMSTEP is omitted, the same check is performed for all simulation steps.
%
%   -3: Running (disabled)
%   -2: Disabled
%   -1: Error in last run (enabled)
%    0: Ready to run
%    1: Completed without errors
%    2: Completed with warnings
%
% See also: SPLITGUI, SETFLAG

    isenabled = @(x) x > -2;
    iscomplete = @(x) x > 0;

    global GUI
    
    allsteps = setdiff(fieldnames(GUI),{'menu'},'stable');
    if nargin < 1 || isempty(simstep), simstep = allsteps; end
    if ~iscell(simstep),simstep = {simstep}; end

    flags = cellfun(@(s) GUI.(s).flag,allsteps);  % vector of all step-flags
    names = cellfun(@(s) get(GUI.(s).btn,'String'),allsteps,'unif',0); % sim-step names
    
    % Recover from previous crash (e.g. onCleanup object calls)
    if any(flags == -3)
        for j = find(flags == -3)'
            setflag(allsteps{j},-1,'Something went wrong, see details in log.');
            flags(j) = -1;
        end
    end
    
    % Update menu
    if ~isempty(getSimOption('prjname'))
        set(GUI.menu.save,'Enable','on'); 
        set(GUI.menu.saveas,'Enable','on');
    end
    if any(GUI.setup.flag == [-1,1,2]), set(GUI.menu.simopt,'Enable','on'); end
    
    if GUI.arrdef.flag > 0 || evalin('base','exist(''ArrayDef'',''var'')'), set(GUI.menu.reduce,'Enable','on'); end
    if GUI.shading.flag > 0 || evalin('base','exist(''ShRes'',''var'')'), set(GUI.menu.export_shdres,'Enable','on'); end
    if GUI.irrtrans.flag > 0 || evalin('base','exist(''EB'',''var'')'), set(GUI.menu.export_balance,'Enable','on'); end
    if GUI.solver.flag > 0 || evalin('base','exist(''SolRes'',''var'')'), set(GUI.menu.export_solres,'Enable','on'); end
    
    if GUI.solver.flag > 0 || GUI.irrtrans.flag > 0, set(GUI.menu.print_summary,'Enable','on'); end
    
    % PROVISIONAL: when all non-computationally-intensive steps are complete...
    if all(cellfun(@(s) GUI.(s).flag,setdiff(allsteps,{'shading','irrtrans','solver'})))
        set(GUI.menu.minprj,'Enable','on'); % enable save-minimal
        % set(GUI.menu.divide,'Enable','on');
        % if ~isempty(dir('.sub/')), set(GUI.menu.merge,'Enable','on'); end
    end

    % Set container-menu Enable status based on Children
    for f = fieldnames(GUI.menu)'
        if isempty(GUI.menu.(f{1}).Children), continue; end
        if any(strcmp('on',{GUI.menu.(f{1}).Children.Enable})), set(GUI.menu.(f{1}),'Enable','on'); 
        else, set(GUI.menu.(f{1}),'Enable','off'); 
        end
    end
        
    % Enable steps when all their requirements are met
    for j = 1:numel(simstep)
        reqok = iscomplete(flags(GUI.(simstep{j}).req));  % bool-vector, required steps completed?
        msg = get(GUI.(simstep{j}).txt,'String');         % current msg in textbox
        if ~iscell(msg), msg = {msg}; end

        if ~isenabled(GUI.(simstep{j}).flag)
            if isempty(reqok)
            % what to do for steps with no requirements... don't do anything for now
            elseif all(reqok)
            % if all requirements are met, set step to ready
                switch simstep{j} % set case-specific message
                case {'shading','irrtrans','solver'}
                    msg = ['Ready to run' names{j}];
                otherwise
                    msg = msg{1};
                end
                setflag(simstep{j},0,msg);
            else
            % otherwise list requirements on second txt line
                missing = names(GUI.(simstep{j}).req);
                missing = missing(~reqok);
                msg{2} = ['Requires ' shortliststr(missing,'',2)];
                setflag(simstep{j},-2,msg); 
                % NOTICE that flag is set to -2 (disabled), even if it was -3 (running)
                % There is no reason why this should happen during operation.
            end
        else
            % Reset all flags, in case things don't match - e.g. when loading an existing project
            setflag(simstep{j},GUI.(simstep{j}).flag,msg);
        end
    end
end
