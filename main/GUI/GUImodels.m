function GUImodels(~,~)
% GUIMODELS - Load (and test) electrical models for modules, diodes, inverter, and cell-temperature

    global GUI
    
    setflag('models',-3,'Runing GUImodels...');
    finalwillandtestament = onCleanup(@updateflags);
    fprintf('\nRuning GUImodels...\n');

    % Load an existing Module-Interpolant, or create one from One-Diode-Model (ODM) parameters
    [ModIVint,ODM] = getmodulemodel();

     % Read Cell-Temperature parameters from ODM and/or SimOptions, create CellTemp interpolant
    [celltemp,NOCT] = getcelltemperaturemodel(ODM,ModIVint);

    % Print module IV curves and basic info
    figh = GUIfigure('modIV','Module Model'); clf();
    ModIVint.plot(celltemp);
    if getSimOption('exportplots'), exportfigure('module',figh); end

    % Generate reverse bias / bypass-diode interpolating structures
    Diode = getbypassdiodemodel(NOCT,ModIVint.Imp0);
    
    % Import Inverter Model, and print
    Inverter = getInverterModel();
    
    [~,~,IAM] = checkIAM([],'-plot');

    % TODO: MODIVINT, IAM, celltemp and Diode should probably be incorporated into a MODULE
    % superclass... maybe once ODM is running on mexfiles and we drop the MODULEINTERPOLANT class.
    
    assignin('base','ModIVint',ModIVint);
    assignin('base','celltemp',celltemp);
    assignin('base','Diode',Diode);
    assignin('base','Inverter',Inverter);
    assignin('base','IAM',IAM);

    % If System is fully defined, check models against it
    ok = true;
    if GUI.arrdef.flag > 0
        ok = checkdesign('-verbose','-soft');
    else
%         warning('GUImodels:notchecked',['Models have not been checked with an array-definition, '...
%             're-run this script just before the Circuit Solution to avoid wasting your time']);
    end
    
    if ok, setflag('models',1); else, setflag('models',2); end
end

function f = findsamlib(libtype)
% Search the current directory for *.samlib files, and return a list of those for which
% LibraryType == libtype

    % Scan all .samlib files for matching libtype
    candidates = cellfun(@dir,{'*.samlib','*.odm'},'unif',false);
    candidates = cat(1,candidates{:});
    seemok = false(size(candidates));
    for j = 1:numel(candidates)
        f = which(candidates(j).name);
        try S = pvlmod_SAMLibraryReader(f); catch, S = struct([]); end
        seemok(j) = isfield(S,'LibraryType') && strcmpi(S.LibraryType,libtype);
    end
    f = {candidates(seemok).name};
end

function [ModIVint, ODM] = getmodulemodel()
% Try to find a saved ModuleInterpolant object

    global GUI
    
    modfile = getmodulefile();
    if strcmpi(modfile(end-5:end),'modint')
        [ModIVint,ODM] = loadexisting(modfile);
    else
        ODM = loadODM(modfile);
        [ModIVint] = createfromODM(ODM);
    end
    if GUI.physmod.flag > 0
        ModIVint = addmodulegeometry(ModIVint);
    end
    
    fprintf(evalc('disp(ODM)'));
    [Pmp,Vmp] = mpp(ModIVint.getIVpp(1000,25));
    fprintf('\tMPP @ STC (1000 W, 25°C): %0.1f W, %0.2f V, %0.2f A\n\n',Pmp,Vmp,Pmp/Vmp);
    
    function modfile = getmodulefile()
        if numel(dir('*.modint'))==1
            modfile = pickfile('*.modint');
        elseif numel(dir('*.modint')) > 1
            modfile = pickfile({'*.modint';'*.samlib','*.odm'},'Select Module Definition file');
        else
            modfile = findsamlib('ODM');
            if numel(modfile)==1, modfile = modfile{1};
            else
                modfile = pickfile({'*.modint;*.samlib';'*.odm'},'Select Module Definition file');
            end
        end
        assert(~isempty(modfile),'GUImodels:nomodfile','No modint/samlib file selected');
    end

    function [ModIVint,ODM] = loadexisting(modfile)
        msg = sprintf('Loading existing ModuleInterpolant from: %s...',relativepath(modfile));
        fprintf('%s\n',msg);
        setflag('models',-3,msg,'app');
        try
            ModIVint = []; %#ok<NASGU>
            load(modfile,'-mat','ModIVint');
            if ~exist('ModIVint','var')|| isempty(ModIVint)||~isa(ModIVint,'ModuleInterpolant')
                clear ModIVint
            else
               % Retrieve ODM parameters, in case the ModuleInterpolant comes from a ODM
               if isstruct(ModIVint.source), ODM = ModIVint.source; 
               else, ODM = struct([]); end
               msg = [ModIVint.name '-Interpolant loaded'];
               fprintf('\t%s\n',msg);
               setflag('models',-3,msg,'cat');
            end
        catch
            error('Could not load selected *.modint file');
        end
    end

    function ODM = loadODM(modfile)
    % Create a new module interpolant, if there isn't one
    
        % Field rename rules for backwards compatibility
        RULES = {'([a-z])0$','$1_ref'; % X0 -> X_ref
                 'Rshbase','Rsh_0';
                 'Rshexp','Rsh_exp'};
    
        msg = sprintf('Reading ODM parameters from : %s...',relativepath(modfile));
        fprintf('%s\n',msg);
        setflag('models',-3,msg,'app');

        ODM = pvlmod_SAMLibraryReader(modfile);
        
        fields = fieldnames(ODM);
        for j = 1:size(RULES,1)
            fields = regexprep(fields,RULES{j,:},'ignorecase');
        end
        if ~isequal(fieldnames(ODM),fields)
            ODM = cell2struct(struct2cell(ODM),fields);
        end

        msg = sprintf('Loaded ODM: %s',ODM.name);
        fprintf('\t%s\n',msg);
        setflag('models',-3,msg,'cat');
    end

    function ModIVint = createfromODM(ODM)
        setflag('models',-3,'Generating Module-Interpolant...','app');

        % Create ModuleInterpolant using getSimOptions & Defaults
        ModIVint = ModuleInterpolant(ODM);
        save([ModIVint.name '.modint'],'ModIVint');

        msg = {[ODM.name ' - Interpolant created'];
               'Performing Interpolation Tests... this can take a few minutes.'};
        setflag('models',-3,msg,-1);
        fprintf('\t%s\n',msg{:});

        InterpolationTest(ModIVint);

        reply = optquestdlg(['Please check the Interpolation Test Results and assert that they''re '...
                'within the required precision. Try using a finer grid (Gegrid, Tcgrid) and/or '...
                'a smaller tolerance (RelTol) to improve accuracy... bear in mind that both of '...
                'these changes will slow down all calculations.'],'Interpolation Test',...
                'They look allright','Quit to fix','They look allright');
        assert(isequal(reply,'They look allright'),'Stopped by user, to fix non-standard Mod. Interp. parameters');

        setflag('models',-3,{'Interpolation Tests complete.'},-1);
        fprintf('\tInterpolation Tests complete.\n');
    end

    function ModIVint = addmodulegeometry(ModIVint)
     % If there is a Module-Interpolant, offer to save the new definition for next time
        try
            ModuleAreas = evalin('base','ModuleAreas');
        catch
            ModuleAreas = evalin('base','Trackers.module.geom');
        end
        if isequal(ModIVint.geom,ModuleAreas), return; end
        
        msg = ['Do you want to save the existing physical-module definition '...
               'into the Module-Interpolant file'];
               
        switch optquestdlg(msg,'GUIphysmod','Yes','No','Yes')
            case 'Yes'
                ModIVint.geom = ModuleAreas;
                modintfile = [strrep(ModIVint.name,'.','_') '.modint'];
                save(modintfile,'ModIVint');
                fprintf('\tPhysical module definition saved to %s.\n', modintfile);
            otherwise
                fprintf('\tPhysical module definition not included.\n');
        end
    end
end

function [celltemp,NOCT] = getcelltemperaturemodel(ODM,ModIVint)
% Read Cell-Temperature parameters from ODM and/or SimOptions, create CellTemp interpolant
    [celltemp,ctf0,S] = getcelltempfcn(ODM,'eff',@ModIVint.MPPefficiency);
    fprintf('Cell-Temp. model: %s\n',S.info);
    setflag('models',-3,['CellTemp: ' S.info],'app');

    if ODM.isCPV, Gnoc = 1000; VWnoc = 2; Tnoc = 20;    % CPV modules, FUTURE: update to standard!              
    else, Gnoc = 800; VWnoc = 1; Tnoc = 20;             % non-concentrating  
    end
    NOCT = ctf0(Gnoc,Tnoc,VWnoc);
    fprintf('\tNOCT: %0.1f °C @ %d W, %d°C, %0.1f m/s\n',NOCT,Gnoc,Tnoc,VWnoc);

    [Pmp,Vmp] = mpp(ModIVint.getIVpp(Gnoc,NOCT));
    fprintf('\tMPP @ NOC: %0.1f W, %0.1f V, %0.2f A\n\n',Pmp,Vmp,Pmp/Vmp);
end

function Diode = getbypassdiodemodel(NOCT,I_ref)
% Generate reverse bias / bypass-diode interpolating structures
    fprintf('Generating Bypass-Diode Interpolant\n');
    Diode = BypassInterpolant(getSimOption('diode')); % use default tolerance and temperatures

    % Plot diode curves, and print basic info in log
    plot(Diode);
    disp(Diode.Pars);
    Vf = inval(Diode.getIVpp(NOCT),I_ref);
    fprintf('\tVf @ %0.1f A, %0.1f °C: %0.3f V\n\n',I_ref,NOCT,Vf);
end

function Inverter = getInverterModel()
% Import Inverter Model

    filename = cat(2,findsamlib('SandiaInverter'),pickfile({'*.ondcp','*.OND'},Inf));
    if numel(filename) ~= 1
        filename = pickfile('*.ondcp;*.samlib;*.OND','Select Inverter Definition file','ui',1);
    else
        filename = filename{1};
    end
    assert(~isempty(filename),'GUImodels:noinvfile','No OND/samlib file selected');
    if ~isempty(regexpi(filename,'.*\.samlib$'))
        fprintf('Loading SANDIA inverter model from *.samlib file...\n');
        Inverter = pvlmod_SAMLibraryReader(filename);
    else
        fprintf('Loading PV-syst inverter model from *.ondcp file...\n');
        Inverter = ONDread(filename);
        Inverter.name = strrep(filename,'\','/'); % \ causes all kinds of scape-char issues
    end

    fh = inverterplot(Inverter);       
    if getSimOption('exportplots'), exportfigure('inverter',fh); end

    fprintf(evalc('disp(Inverter)'));
    fprintf('\b');
    setflag('models',-3,['Inverter: ' Inverter.name],'app');
end
