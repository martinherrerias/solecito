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
    
    % Load custom options from PRJNAME_simoptions.json, complete with defaults
    optionsfile = regexprep(SimOptions.prjname,'(\.ssp)$','*.json');
    optionsfile = pickfile(optionsfile,1,'options.json file (cancel to create default)','-soft');
    if ~isempty(optionsfile)
        try
            readoptionsfile(optionsfile);
            flagmsg{1} = sprintf('Using %s',optionsfile);
            fprintf('\tUsing %s...\n',optionsfile);
        catch ERR
            warning('Failed to read %s, custom options ignored! - error: %s\n',...
                optionsfile,getReport(ERR));
            optionsfile = {};
        end
    end
    if isempty(optionsfile)
        optionsfile = writeoptionsfile();
        fprintf(['\t%s options file created in project directory, ',...
                 'make changes and run Setup again, if required.\n'],optionsfile);
        flagmsg{1} = [optionsfile ' not found, creating default'];
    end

    allflags = nestedstruct2cell(OptionFlags);
    allflags = [allflags{:}]';
    flagmsg{2} = sprintf('%d of %d option fields match DefaultOptions', nnz(allflags), numel(allflags));

    evalin('base','global SimOptions;'); 
    evalin('base','global OptionFlags;');

    % assignin('base','diode',SimOptions.diode);
    % assignin('base','CellTemp',SimOptions.CellTemp);
        
    printsummary('setup');
    setflag('setup',1,flagmsg);
end
