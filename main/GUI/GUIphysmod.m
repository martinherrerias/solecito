function GUIphysmod(varargin)
% Read physical module definition - see stdmodulepoly() for file generation

    global GUI

    setflag('physmod',-3,'Runing GUIphysmod...');
    finalwillandtestament = onCleanup(@updateflags);
    fprintf('\nRunning GUIphysmod...\n');
    
    Trackers = evalin('base','Trackers');

    if nargin == 1
    % Allow reuse of existing Mod = Trackers.module
        Mod = varargin{1};
        gotitfrom = 'prev';
        Trackers.module = Mod;
    else
        if isstruct(Trackers.module), Trackers.module = Trackers.module.default; end

        % Set module-definition source priority: modID.mpoly > modID.modint > user-pick
        % Use of WHICH(..) allows having a module-library included in PATH
        filename = which([Trackers.module '.mpoly']);
        if isempty(filename), filename = which([Trackers.module '.modint']); end
        if isempty(filename)
        % Aks the user to pick a (possibly non-matching) *.mpoly / *.modint file
            msg = sprintf('Cannot find module-geometry-file: %s.mpoly. ',Trackers.module);
            msg = [msg,'Do you want to search for it elsewhere (optionally selecting a ',...
                       'replacement-model), or attempt to generate a standard-cell-size ',...
                       'module that fits the provided mount?'];
            switch questdlg(msg,'Missing module-file','Browse/Replace','Standard','Standard')
            case 'Browse/Replace'
                filename = pickfile({'*.mpoly';'*.modint'},'Select module-geometry-file','ui',1);
            case 'Standard'
            % Try to create a standard-module based on Trackers.geom.elements
                filename = [Trackers.module '.mpoly'];
                stdmodulepoly(Trackers.geom.elements(1),'outputfile',filename);
                setflag('physmod',-3,'Creating standard-cell physical module...');
                fprintf('Generating standard-cell-size module from mount-element dimensions.\n');
                fprintf('\tIf you don''t like what you get, create your own *.mpoly files and re-run this script.\n');
                fprintf('\tUse stdmodulepoly() to generate custom rectangular module files.\n');
                gotitfrom = 'std'; % standard
            otherwise
               filename = 0;
            end
        end

        % Attempt to actually read whatever source was selected
        if isequal(filename,0)
            Mod = [];
        else
            setflag('physmod',-3,'Reading physical module definition...');
            fprintf('Reading physical module definition from: %s\n',filename);
            if numel(filename) > 7 && strcmp(filename(end-6:end),'modint')
            % Get module-geometry from module-interpolant (*.modint)
                Mod = load(filename,'-mat','ModIVint'); 
                Mod = Mod.ModIVint;
                gotitfrom = 'interp'; % interpolant
            else
            % Get module-geometry from model-geometry file (*.mpoly)
                Mod = importpolygonsfile(filename);
                gotitfrom = 'mpoly'; % file
            end
        end
        if ~isfield(Mod,'geom'), Mod = []; end
        assert(~isempty(Mod),'GUIphysmod:nomod','Failed to create module-geometry');

        Trackers.module = struct('default',Trackers.module,'loaded',Mod.name,'source',filename,'geom',Mod.geom);
        if isfield(Mod,'info'), Trackers.module.info = Mod.info; end
    end

    % Attempt to replicate modules into mount, issue warnings/error if they don't fit
    Trackers.geom = fitmoduleinmount(Mod.geom,Trackers.geom);

    assignin('base','Trackers',Trackers);

    switch gotitfrom
        case 'mpoly', flagmsg{1} = 'Physical module imported.';
        case 'interp', flagmsg{1} = 'Physical module read from Module-Interpolant.';
        case 'std', flagmsg{1} = 'Standard Physical module created.';
        case 'prev', flagmsg{1} = 'Physical module preserved.';
    end

    Nc = prod(Mod.geom.dims);
    Nd = Mod.geom.dims(1);
    [w,h,~] = rectangleproperties(Mod.geom.border);
    flagmsg{2} = sprintf('%d cell, %d cell-string, (%0.1f x %0.1f cm) module',Nc,Nd,w*100,h*100);
    fprintf('\t%s\n',flagmsg{2});
    setflag('physmod',1,flagmsg);

    assignin('base','ModuleAreas',Mod.geom);

    % If there is a Module-Interpolant, offer to save the new definition for next time
    if ~strcmp(gotitfrom,'interp') && GUI.models.flag > 0 && isempty(evalin('base','ModIVint.geom'))
        ModIVint = evalin('base','ModIVint');
        msg = sprintf(['Do you want to save the physical-module definition you just '...
                       'created into the current Module-Interpolant?']);
        switch optquestdlg(msg,'GUIphysmod','Yes','No','No')
            case 'Yes'
                ModIVint.geom = Mod.geom;
                modintfile = [ModIVint.name,'.modint'];
                save(modintfile,'ModIVint');
                fprintf('\tPhysical module definition saved to %s.\n', modintfile);
            otherwise
                fprintf('\tModule-Interpolant not modified.\n');
        end
    end

    figh = GUIfigure('trckareas','Physical Module'); clf(figh)
    labeling = prod(Trackers.geom.dims(1:2)) < 100;
    plotArea(Trackers.geom,3,labeling);
    plotArea(Trackers.geom.elements(1));
    set(gca,'Color',[0.8;0.8;0.8]); 
    axis equal; axis tight; title('Detailed Mount'); xlabel('x [m]'); ylabel('y [m]');
    hold on
    plot(Trackers.analysedpoints(1,:),Trackers.analysedpoints(2,:),'kx','MarkerSize',10)

    if getSimOption('exportplots'), exportfigure('physmod',figh); end
    
end

