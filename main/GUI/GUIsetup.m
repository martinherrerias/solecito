function GUIsetup(~,~)
    % PROVISIONAL: should be replaced by one or more configuration files
    % Structure fields in: SimOptions, diode, CellTemp, Tracker
    % PENDING: Clarify separation between Default and Custom options
    
    global SimOptions;
    global OptionFlags; % copy of SimOptions with boolean fields 
                        % OptionFlags.(j) is true if SimOptions.(j) wasn't set explicitly
                        
    setflag('setup',-3,'Runing GUIsetup...');
    fprintf('\nSimulation setup:\n');
    finalwillandtestament = onCleanup(@updateflags);
            
    % Clear persistent variables in other functions (used to avoid repetitive getSimOption calls)
    clearvars -global -except GUI SimOptions OptionFlags
    clear PWLparallel PWLseries mkivpp parsetolerance % clear persistent
    clear stopwatch radiance_igawa onediodeMPP onediodemodel onediodemodel2
    
    % Clear SimOptions structure
    fields2keep = {'version','prjname'};
    SimOptions = rmfield(SimOptions,setdiff(fieldnames(SimOptions),fields2keep));
    OptionFlags = cell2struct(repmat({true},[numel(fields2keep),1]),fields2keep);
    
    % First guess is PRJPATH/AddSimOptions.m
    optionsfile = fullfile(fileparts(SimOptions.prjname),'AddSimOptions.m');
    if isempty(dir(optionsfile))
        if ~isempty(dir('./AddSimOptions.m'))
            copyfile('./AddSimOptions.m',optionsfile);
            fprintf(['\tPWD/AddSimOptions.m file copied to project directory, ',...
                 'make changes and run Setup again, if required.\n']);
            flagmsg{1} = 'AddSimOptions.m copied from Working directory';
        else
            writeoptionsfile(fileparts(optionsfile));
            fprintf(['\tDefault AddSimOptions.m file added to project directory, ',...
                     'make changes and run Setup again, if required.\n']);
            flagmsg{1} = 'AddSimOptions.m not found, using defaults';
        end
    else
        fprintf('\tUsing PRJPATH/AddSimOptions.m...\n');
        flagmsg{1} = 'Using PRJPATH/AddSimOptions.m';
    end
    
    fprintf('\tExecuting AddSimOptions.m...\n');
    run(optionsfile);
    
    fprintf('\tCompleting SimOptions structure...\n');
    completeoptions(); % update/complete SimOptions and OptionFlags with DEFAULTOPTIONS

    allflags = nestedstruct2cell(OptionFlags);
    allflags = [allflags{:}]';
    flagmsg{2} = sprintf('%d of %d option fields completed from DefaultOptions', nnz(allflags), numel(allflags));

    evalin('base','global SimOptions;'); 
    evalin('base','global OptionFlags;');

    % assignin('base','diode',SimOptions.diode);
    % assignin('base','CellTemp',SimOptions.CellTemp);
        
    printsummary('setup');
    setflag('setup',1,flagmsg);
end
