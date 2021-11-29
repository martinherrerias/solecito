function SolRes = pvArraySolver(varargin)
% SolRes = pvArraySolver(POA,ShRes,ArrDef,ModIV,Diode,Inverter,[OPT,'opt',val,..])
%   Calculate PV-plant's DC operating points over a time-series, given irradiance and temperatures 
%   for each tracker at every step (POA), electrical models for all system's components, a connection-
%   scheme for the complete plant, and 
%   beam-shading information for every tracker and time-step.
%
% SolRes = pvArraySolver('backup',S) - attempt to resume interrupted calculation from backup
%   structure S = load(FILE), after a previous call with ..,'backup',FILE.
%
% INPUT:
%   POA (Plane-Of-Array) structure containing irradiance and cell-temperature information for 
%       every mount (1:Nu) and module (1:Nm) at every simulation time-step (t = 1:Nt). 
%       See EFFECTIVEIRRADIANCE(...) for details:
% 	POA.Dpoa - Nt·Nu·Nm array with intensity of diffuse irradiance components after shading, 
%       soiling, spectral-correction, Incidence-Angle-Modified-Absortivity, and intra-module
%       diffuse mismatch for every MODULE in the array.
% 	POA.Bpoa - Nt·Nu* array with intensity of direct irradiance after diffuse-shading, soiling, 
%       spectral-correction, and Incidence-Angle-Modified-Absortivity. SHADING NOT INCLUDED
%   POA.ShF - Nt·Nu·Nm array with fractions of shade-coverage for each module, mount, and time-
%       step. ("Linear" shading, 1 = full shade, 0 = no shade)
%   POA.Tush - Nt·Nu·Nm array of cell-temperatures for modules which are NOT shaded
%   POA.Tsh - Nt·Nu·Nm array of cell-temperatures for shaded modules
%   POA.Nsb - Nt·Nu·Nm array, number of shaded cell-blocks on every module, mount, and timestep.
%   POA.Defects - (optional) [Nu·Np,Nc] constant, cell-specific irradiance multiplying factors
%       to simulate manufacturing & degradation mismatch. See RANDOMDEFECTS.
%
%   ShRes (SHADINGRESULTS object) contains (among others) the results of a Beam-Shading analysis, namenly:
%   ShRes.BLPidx - Nt·Nat array of indices for beam-light polygons. Serves only to avoid storing empty polygons for all
%       trivial cases (i.e. no shade, and full-shade):
% 			BLPidx(t,tr) = 0 means there is no shading over tracker tr at timestep t
% 			BLPidx(t,tr) = -1 means tracker tr is fully shaded at timestep t
% 			BLPidx(t,tr) = k, k > 0 implies partial shading, with beam-sunlight over polygon ShRes.BLpoly{k}
%   ShRes.BLpoly - Cell-Array of beam-light polygons (i.e. complement of beam shades over the tracker plane)
%       Each element can be a single object, or an array of objects of the polygon class.
%
%   ShRes.MountGeom - a pvArea object representing a tracker and its subdivisions:
%       MountGeom.border - a polygon object representing the tracker outline
%       MountGeom.elements(p) - a pvArea object representing module p in the tracker
%       MountGeom.elements(p).elements(b) - cell-block b inside module p
%       MountGeom.elements(p).elements(b).elements(c) - cell c inside cell-block b
%     Each subdivision contains also an outline polygon (in tracker-centered-coordinates) such that 
%     e.g. MountGeom.elements(p).elements(b).border returns the outline of cell-block b.
%
%   ArrDef - NDmap object which includes the connection scheme for all modules in the plant, with
%       E-coordinates [MPPT-input, string-number, position-in-string] and P-coordinates [mount,
%       module-position].
%
%   ModIV - ModuleInterpolant class object. Defines a PV-module's behavior from a set of reference curves at diferent 
%       effective irradiance and cell temperature conditions.
%
%   Diode - BypassInterpolant class object. Defines a reverse-biassed Diode's behaviour from a set of ideal parameters.
%
%   Inverter - Parameter structure defining the electrical model of an inverter. Used only to define expected limits for
%       array voltage and current (might be changed in future versions).
%
% OPT - options structure, provided as such (and completed with default SimOptions), and/or
%   corrected by individual 'name',value pairs:
%
%   OPT.bypassdiodes - bypass diodes in every cell-block
%   OPT.blockingdiodes - diodes in series with every string
%   OPT.parallelelements - cell-blocks are connected in parallel, not series
%   OPT.modulediodes - bypass diode for complete module
%
%   OPT.cellshadingtol - smallest difference to consider between two shaded cells
% 	OPT.quickshading - Approx. factors based on the number of shaded cell-border vertices
%   OPT.skipnonshaded - override Gvartol, calculate mismatch for all PV elements (slow!)
%                           consider reducing Gvartol before setting this to false
%   OPT.celldefects - consider or ignore POA.Defects
%
% 	OPT.Gvartol - Irradiance variability within a string above which circuit must be solved
% 	OPT.RelTol - Tolerance for PWL curve simplification
%   OPT.Lims_Voc0 - Limits for IV curve (voltage in Voc0 units)
%   OPT.Lims_Isc0 - Limits for IV curve (voltage in Isc0 units
%   OPT.MPPmethod - {'interp','ivpp','odm'} used by ModuleInterpolant.getMPP for uniform cases
%
%   OPT.exportplots - export each plot to a PNG file
%   OPT.solverlog - use dedicated (verbose) log
% 	OPT.plotting - plot array IV curve for shaded trackers (slow!)
% 	OPT.savingIVpps - save all array IVpp's (memory!)
%
%   OPT.backup - Name for auto-save dump -mat file
%   OPT.autosavemin - auto-save workspace every x minutes
%
% OUTPUT:
% 	SolRes.Push - Nt·Ni·Ns array. Push(t,i,s) returns the added MPP power of all modules of string s on MPPT i at
%       time-step t, under operating conditions, but as if disconnected from other strings, and not-shaded.
% 	SolRes.Vush - Nt·Ni·Ns array. MPP string voltages, for the conditions above.
% 	SolRes.Psh - Nt·Ni·Ns array. MPP of every string, considering beam shading over its modules, and missmatch due 
%       to connection with other strings of the MPPT-input.
% 	SolRes.Vsh - Nt·Ni array. MPP voltage for the above conditions (all strings in an MPPT share a common voltage)
% 	SolRes.Pdc - Nt·Ni·Ns array. MPP of every string, considering beam shading over its modules, missmatch due 
%       to connection with other strings, and Inverter Limitations (clipping).
% 	SolRes.Vdc - Nt·Ni array. MPP voltage for the above conditions (all strings in an MPPT share a common voltage).
% 	SolRes.simtimer - stopwatch class object with simulation-timing information.

%% Parse input, and complete simulation options with defaults

    global SimOptions
    SimOptions = completeoptions(); % fill with defaults if necessary

    % Required simulation options
    REQ_OPT = {'Gvartol','Lims_Isc0','Lims_Voc0','MPPmethod','RelTol','autosavemin',...
        'blockingdiodes','bypassdiodes','cellshadingtol','exportplots','modulediodes',...
        'parallelelements','plotting','quickshading','savingIVpps','skipnonshaded','solverlog',...
        'backup','celldefects'};
    opt = rmfield(SimOptions,setdiff(fieldnames(SimOptions),REQ_OPT));
    clear SimOptions
    opt.backup = '';
    [opt,varargin] = getpairedoptions(varargin,opt);

    % Mininum list of variables that allow resuming after crash
    VARS_TO_SAVE = {'POA','ShRes','ArrDef','ModIV','Diode','Inverter','opt',...
                    'Push','Vush','Psh','Vsh','Pdc','Vdc','Plsh','Pwcs','simtimer','t','savedIVpps'};

    switch numel(varargin) 
    case 0
    % SolRes = pvArraySolver('backup',S)
    % Resume from crash, load existing results from previous partial simulation structure...
        assert(isstruct(options.backup) && all(isfield(options.backup,VARS_TO_SAVE)),...
            'Failed to recover from backup');

        varargin = struct2cell(orderfields(options.backup,VARS_TO_SAVE));
        [POA,ShRes,ArrDef,ModIV,Diode,Inverter,opt,...
            Push,Vush,Psh,Vsh,Pdc,Vdc,Plsh,Pwcs,simtimer,t0] = deal(varargin{:});

        t0 = t0 + 1; % start at the next not-finished timestep
    case {6,7}
    % SolRes = pvArraySolver(POA,ShRes,ArrDef,ModIV,Diode,Inverter,[OPT])
        if numel(varargin) == 7
            opt = completestruct(varargin{end},opt);
            opt = rmfield(opt,setdiff(fieldnames(opt),REQ_OPT));
        end
        [POA,ShRes,ArrDef,ModIV,Diode,Inverter] = deal(varargin{1:6});
        
        simtimer = stopwatch();
        t0 = 1;
    end
    clear varargin

    % Parse POA structure, expand Bpoa to be [Nt,Nu,Np]
    parsestruct(POA,{'Dpoa','Bpoa','Tush','Tsh','ShF','Nsb'},'-r','-n');
    [POA.Dpoa,POA.Bpoa] = compatiblesize(POA.Dpoa,POA.Bpoa,POA.Tush,POA.Tsh,POA.ShF,POA.Nsb);
    POA = rmfield(POA,setdiff(fieldnames(POA),{'Dpoa','Bpoa','Tush','Tsh','ShF','Nsb'}));
    
    POA.Gpoa = POA.Dpoa + (1-POA.ShF).*POA.Bpoa;            % Linearly shaded!
    POA.Tavg = POA.Tsh.*POA.ShF + (1-POA.ShF).*POA.Tush;    % id.

    notok = POA.Dpoa + POA.Bpoa > ModIV.Gref(end);
    if any(notok,'all')
        warning('%d module GTI values at %d time-steps > %0.1f W/m² (interpolant limit) set to NaN',...
            nnz(notok),nnz(any(notok,2:3)),ModIV.Gref(end));
        POA.Dpoa(notok) = NaN;
        POA.Bpoa(notok) = NaN;
        POA.Gpoa(notok) = NaN;
    end
    notok = POA.Dpoa < 0 | POA.Bpoa < 0 | POA.Gpoa < 0;
    if any(notok,'all')
        warning('%d negative module GTI values at %d time-steps set to NaN',nnz(notok),nnz(any(notok,2:3)));
        POA.Dpoa(notok) = NaN;
        POA.Bpoa(notok) = NaN;
        POA.Gpoa(notok) = NaN;
    end
    notok = POA.Tush < ModIV.Tref(1) | POA.Tush > ModIV.Tref(end) | ...
            POA.Tsh < ModIV.Tref(1) | POA.Tsh > ModIV.Tref(end);
    if any(notok,'all')
        warning('%d cell-temperature values at %d time-steps set to NaN',nnz(notok),nnz(any(notok,2:3)));
        POA.Dpoa(notok) = NaN;
        POA.Bpoa(notok) = NaN;
        POA.Gpoa(notok) = NaN;
    end
    [Nt,Nu,Np] = size(POA.Dpoa); % No. of Time-Steps, No. of mounts, No. modules per mount

    % Parse shading-results
    if isa(ShRes,'ShadingResults')
    % ... keep only beam-shading-related indices & polygons
        ShRes = struct(ShRes,'except',{'info','worldgeom','BshF','Nsb','DwF0','DshF','DwF'});
    end
    parsestruct(ShRes,{'fullshaded','partshaded'},'-l','size',[Nt,Nu]);
    parsestruct(ShRes,{'partidx'},'-i','size',[Nt,Nu]);
    parsestruct(ShRes,{'BLPidx','BLPoly','mountgeom'});
    assert(issparse(ShRes.BLPidx) && size(ShRes.BLPidx,1) == nnz(ShRes.partshaded) && ...
        size(ShRes.BLPidx,2) == numel(ShRes.BLPoly),'Bad polyon indices/shading-polygons');

    % ... parse mount-geometry, while initializing cell-shading-calculation function
    [~,Nepm,Ncpe] = CellShadingFactors(ShRes.mountgeom,Np);

    % Parse array-definition
    assert(isa(ArrDef,'NDmap') && ArrDef.psize(1) == Nu && ArrDef.psize(2) == Np,...
        'Expecting [%d %d]-psize NDmap ArrDef',Nu,Np);

    Ni = ArrDef.esize(1); % No. of inverters (MPPT-inputs)
    Ns = ArrDef.esize(2); % Max. No. of strings per MPPT 
    Nm = ArrDef.esize(3); % Max. No. of modules per string

    % Parse Module- and Diode-Interpolants
    assert(isa(ModIV,'ModuleInterpolant'),'');
    assert(isa(Diode,'BypassInterpolant'),'');

    % Parse inverter model
    parsestruct(Inverter,{'MPPTLow','MPPTHi','Ps0','Pdc0'},'-n','-f','-r','-p');
    Inverter = rmfield(Inverter,setdiff(fieldnames(Inverter),{'MPPTLow','MPPTHi','Ps0','Pdc0','IdcMax'}));

    plotobj = arraysolverplotter(false);  % start with 'waitbar' instance of plotter    
    cleanupobj = onCleanup(@() plotobj.close(true)); % close on error
    if plotobj.isUI, setbtn(plotobj,[],'Enable','off'); end  % disable 'switch-to-plotter' button
    plotobj.updatewaitbar(0,'Parsing input data...');

    % labels are used in log and plots, try to get them from POA/ShRes, but don't sweat it
    try timelabels = datestr(POA.t,'yyyy/mm/dd HH:MM'); 
    catch, timelabels = repmat('?',Nt,1);
    end
    try stepinfo = arrayfun(@(j) sprintf('Az:%0.1f°, El:%0.1f°, avg(Gpoa):%0.0f W/m²',...
                            ShRes.az(j),ShRes.el(j),mean(POA.Gpoa(j,:))), 1:Nt,'Unif', false)';
    catch, stepinfo = repmat({'?'},Nt,1);
    end

    if t0 == 1
    % Allocate result arrays
        Push = zeros(Nt,Ni,Ns);		% DC max. power (sum of un-shaded individual modules)
        Vush = zeros(Nt,Ni,Ns);
        Psh = zeros(Nt,Ni,Ns);		% DC max. power after shading & mismatch
        Vsh = zeros(Nt,Ni);
        Pdc = zeros(Nt,Ni,Ns);		% DC power subject to MPPT bounds
        Vdc = zeros(Nt,Ni);

        Plsh = zeros(Nt,Ni);    % Linear shading
        Pwcs = zeros(Nt,Ni);    % Worst-case shading
        
        if opt.savingIVpps, savedIVpps(Nt,Ni) = mdpwl(); % allocate an array of MDPWL curves
        else, savedIVpps = mdpwl.empty();
        end 
    end

    % Set expected I, V, and P limits
    %arrIVLimits = [-Inf Inverter.Vdcmax -Inf Inverter.Idcmax];
    %Tol = [Inverter.Vdc0, Inverter.Pdc0/Inverter.Vdc0]*opt.RelTol;
    %modIVLimits = [-Inf 1.2*ModIV.Voc0 -Inf 1.5*ModIV.Isc0];
    modIVLimits = [opt.Lims_Voc0*ModIV.Voc0,opt.Lims_Isc0*ModIV.Isc0];
    arrIVLimits = modIVLimits.*[1 Nm 1 Ns];
    Tol = [ModIV.Vmp0, ModIV.Imp0 1]*opt.RelTol; % Goes to addseries/parallel/array -> MDPWL
    Plim = [Inverter.Ps0,Inverter.Pdc0];
    Vlim = [Inverter.MPPTLow,Inverter.MPPTHi];

    simtimer.add2lap('parsing');

    % Find the minimum irradiance value for which it is worth performing the calculations:
    % Get the array's max. production at Gmin and min(Tc), and iterate to find P(Gmin) = Ps0
    Tmin = max(min([POA.Tush(:);POA.Tsh(:)]),ModIV.Tref(1));
    Pthresh = @(g) ModIV.getMPP(g,Tmin,opt.MPPmethod)*Ns*Nm - Plim(1);
    Gmin = bisection(Pthresh,0,500,0.1);
    Gmin = Gmin - 0.1; % ... and make sure we're below
    fprintf('\tIrradiance threshold set at %0.1f W/m²\n',Gmin);

    plotobj.updatewaitbar(0,'Creating cell-interpolant...');
    % Get cell IV interpolant from Module IV interpolant
    if opt.parallelelements
        cellIV = scale(ModIV,1/Ncpe,1/Nepm);
    else
        cellIV = scale(ModIV,1/(Ncpe*Nepm),1);
    end

    % Get the number of bypass-diodes in a module 'as seen from the outside'
    % dpm = [diodes in series, diodes in parallel]
    diodespermodule = [0,0];
    if opt.modulediodes
        if opt.parallelelements || Nepm == 1
        % Nepm + 1 diodes in parallel... not the best design idea, but hey, who am I to judge?
            diodespermodule = [1,Nepm+1];  
        else
        % assume only the 'outer' bypass diode is 'visible' from the outside
            diodespermodule = [1,1];       
        end
    elseif opt.bypassdiodes
        if opt.parallelelements, diodespermodule = [1,Nepm];  % Nepm diodes in parallel
        else, diodespermodule = [Nepm,1];                     % Nepm diodes in series
        end
    end

    % TEST: RANDOM CELL DEFECTS
    if isfield(POA,'Defects')
        assert(all(size(POA.Defects)==[Nu*Np,Nepm*Ncpe]),'POA.Defects size does not match [NmNp,NeNc]');
        % Get [Nu,Np] boolean index of defective modules
        defective = reshape(any(POA.Defects,2),[Nu,Np]);

        % Get [Ni,Ns] boolean index of strings with defective modules
        hasdefective = false(Ni,Ns,Nm);
        hasdefective(ArrDef.eidxT(defective')) = true;
        hasdefective = any(hasdefective,3);
    else
        defective = false(Nu,Np);
        hasdefective = false(Ni,Ns);
    end
    opt.celldefects = any(hasdefective(:));
 
    plotobj.updatewaitbar(0,'Identifying mismatch conditions...');

    usedstrings = false(Ni,Ns);
    dark = true(Nt,Ni,Ns);          % Individual strings below threshold irradiance
    missing = true(Nt,Ni,Ns);       % Individual strings with any NaN values
    mismatched = false(Nt,Ni,Ns);

    mps = shiftdim(sum(ArrDef.PIDX>0,1),1); % an [Ns,Ni] matrix with N° of modules per string
    stringspermppt = sum(mps > 0,1)';       % an [Ni,1] vector
    modulesperstring = max(mps,[],1)';      % an [Ni,1] vector - constant mps for the same mppt

    if any(any(mps > 0 & mps ~= repmat(modulesperstring',Ns,1)))
       error('pvArraySolver:mps','There are strings with different number of modules'); 
    end
    clear mps

    % Create filters to be used later, to avoid unnecesary calculations
    for mppti = 1:Ni
        usedstrings(mppti,:) = ArrDef.PIDX(1,:,mppti) > 0; % look for first module of every string
        for str = find(usedstrings(mppti,:))
            % get physical indices for modules of string k in MPPT-input j
            in_str = ArrDef.pidxT(mppti,str,:);
            
            % Get that string module's irradiances, temperatures, and shading factors
            Dstr = POA.Dpoa(:,in_str);    % diffuse tilted irradiance (shaded)
            Bstr = POA.Bpoa(:,in_str);    % beam tilted irradiance (not shaded)
            ShFj = POA.ShF(:,in_str);     %
            Gstr = Dstr+(1-ShFj).*Bstr;   % Avg. module GTI
            Ustr = Dstr+(ShFj < 1).*Bstr; % GTI on the cell with max. irradiance

            missing(:,mppti,str) = ~all(isfinite(Gstr),2);
            dark(:,mppti,str) = max(Ustr,[],2) < Gmin;

            % Consider as candidates for mismatch modules with |dg/g| > opt.Gvartol
            % Since Gstr is the linearly-shaded Gpoa, this includes also any relevant
            % [i.e. Bpoa > tol] combination of fully-shaded & unshaded modules:
            m = any(abs(1-Gstr./mean(Gstr,2)) > opt.Gvartol,2);
            
            % Include also any strings with partially-beam-shaded modules
            m = m | any(ShFj > 0 & ShFj < 1 & Bstr./Gstr >= opt.Gvartol,2);
            
            if ~opt.skipnonshaded, m = true; end
            mismatched(:,mppti,str) = m & ~dark(:,mppti,str) & ~missing(:,mppti,str);

            plotobj.updatewaitbar((str+Ns*(mppti-1))/(Ni*Ns),'Pre-analysis...');
        end
        
        % Set values for non-existing strings to NaN
        Push(:,mppti,~usedstrings(mppti,:)) = NaN;
        Vush(:,mppti,~usedstrings(mppti,:)) = NaN;
        Psh(:,mppti,~usedstrings(mppti,:)) = NaN;
        Pdc(:,mppti,~usedstrings(mppti,:)) = NaN;
    end
    missingMPPT = any(missing,3);
    hasdefective = hasdefective & usedstrings;
    
    % Flag any defective strings as possibly mismatched
    mismatched = mismatched | (shiftdim(hasdefective,-1) & ~dark);
    progresscount = [nnz(mismatched(1:t0-1,:,:)),nnz(mismatched)];

    % Set missing values to NaN on result structures...
    % ... String-independent [Nt,Ni,Ns] arrays
    Push(missing) = NaN; % DC max. power (sum of un-shaded individual modules)
    Vush(missing) = NaN;

    % ... MPPT-constant [Nt,Ni] arrays
    Vsh(missingMPPT) = NaN;
    Vdc(missingMPPT) = NaN;
    Plsh(missingMPPT) = NaN;    % Linear shading
    Pwcs(missingMPPT) = NaN;    % Worst-case shading

    % ... MPPT-dependent [Nt,Ni,Ns] arrays
    Psh(repmat(missingMPPT,1,1,Ns)) = NaN; % DC max. power after shading & mismatch
    Pdc(repmat(missingMPPT,1,1,Ns)) = NaN; % DC power subject to MPPT bounds

%% Main Loop

    % Prepare full-plotter (just in case)
    
    if plotobj.isUI, plotobj.setaxislimits(polygon.unpack(ShRes.mountgeom.border),Inverter); end
    if opt.plotting
        plotobj.switchtoplotter(); % switch to full-plotter
        plotobj.pausing = true;
    else
        plotobj.fighandle.reset(); 
        if plotobj.isUI, plotobj.setbtn([],'Enable','on'); end  % enable 'switch-to-plotter'
    end

    simtimer.add2lap('preprocessing');
    simtimer.resetlocal();
    lastsaved = simtimer.globaltime();
    pvsLog = NaN; % Used for step-by-step output
    fibo = [0 1]; % DEBUG: save simtimer splits @ t = 1,2,3,5,8,13...

    strIVpp(Ns,1) = mdpwl(); % will hold the IV curves of all strings of an MPPT
    modIVpp(Nm,1) = mdpwl(); % will hold the IV curves of all modules in one string

    for t = t0:Nt

        if t-t0 == sum(fibo)
          simtimer.addsplit(sprintf('main_loop_%d',sum(fibo)));
          fibo = fibo(2) + [0 fibo(1)];
        end

        for mppti = 1:Ni

            % if darkMPPT(t,mppti) || missingMPPT(t,mppti), continue; end % it's 3:00 am, go home

            used = usedstrings(mppti,:);
            if ~any(used), continue; end
            tricky = shiftdim(mismatched(t,mppti,:),1) | hasdefective(mppti,:);

            % Mean MPPT irradiance and temperature
            allmodules = ArrDef.pidxT(mppti,:,:);
            Garr_avg = mean(double(POA.Gpoa(t,allmodules)));
            Tarr_avg = mean(double(POA.Tavg(t,allmodules)));
            
            if ~isfinite(Garr_avg)
            % invalid irradiance data
                Push(t,mppti,used) = NaN;   Vush(t,mppti,used) = NaN;         
                Psh(t,mppti,used) = NaN;    Vsh(t,mppti,1) = NaN; 
                Pdc(t,mppti,used) = NaN;    Vdc(t,mppti) = NaN;
                Plsh(t,mppti) = NaN;
                Pwcs(t,mppti) = NaN;
                continue;
            elseif Garr_avg <= 0
                continue; % let all results in zeros
            end
            
            % Get a single (equal) curve with array's mean irradiance and temperature
            strIVpp(1) = scale(ModIV.getIVpp(Garr_avg,Tarr_avg),modulesperstring(mppti),1);
            [Pavg,Vavg] = strIVpp(1).mpp();
            Plsh(t,mppti) = Pavg*nnz(used);
            
            [Gwc,idx] = min(double(POA.Dpoa(t,allmodules)+(POA.ShF(t,allmodules)==0).*POA.Bpoa(t,allmodules)));
            Twc = POA.Tavg(t,allmodules(idx));
            % Twc = max(double(POA.Tush(t,allmodules)));
            Pwcs(t,mppti) = ModIV.getMPP(Gwc,Twc,opt.MPPmethod)*numel(allmodules);
            
            % if irradiance is uniform (or too low)...
            if ~any(tricky)

                Push(t,mppti,used) = Pavg;
                Vush(t,mppti,used) = Vavg;
                Psh(t,mppti,used) = Pavg;
                Vsh(t,mppti,1) = Vavg; 

                % Find Pmp and Vmp within inverter MPPT limitations
                [Pdc(t,mppti,used),Vdc(t,mppti)] = strIVpp(1).mpp(Vlim);
                continue;
            end
            simtimer.add2lap('array');

            strIVpp(:) = mdpwl();
            for str = find(used)

                % get module indices for string k in MPPT-input j
                % (indices are unique module-numbers, for the complete plant)
                in_str = ArrDef.pidxT(mppti,str,:);
                nm = numel(in_str);
                
                Gsh = double(POA.Dpoa(t,in_str));          % beam-shaded irradiance on str
                Gush = Gsh + double(POA.Bpoa(t,in_str));   % not-beam-shaded ...
                Tush = double(POA.Tush(t,in_str));
                %Tsh = double(POA.Tsh(t,in_str));
                %kd = Gsh./Gush;
                %kd(Gush == 0) = 1;
                
                % Get shading factors for all modules in the current string
                ShFj = double(POA.ShF(t,in_str));

                % Calculate non-shaded, non-clipped MPP
                [Pjk,Vjk] = ModIV.getMPP(Gush,Tush,opt.MPPmethod);
                Push(t,mppti,str) = sum(Pjk);

                % ... assume an ideal buck-converter is used to step-down voltage at every module, 
                % to match the highest MPP current: Is = max(Pjk./Vjk)
                Vush(t,mppti,str) = sum(Pjk)/max(Pjk./Vjk);

                % Now get mean module-average irradiance and temperature (considering shading)
                Gavg = Gush.*(1-ShFj) + Gsh.*ShFj;
                Tavg = double(POA.Tavg(t,in_str));
                
                Gs = mean(Gavg(:));   % Average irradiance for all modules
                Ts = mean(Tavg(:));   % Average cell-temperature ...
                      
                % Add the string to the MPPTi's linear and worst-case power
                % NOTE: linear losses on the beam component only (proportional to 1-kd)
                % Plsh(t,mppti) = Plsh(t,mppti) + sum((1-ShFj).*Pjk.*(1-kd)+kd.*Pjk,2);
                % worstcaseShFj = POA.Nsb(t,in_str) > 0;
                % Pwcs(t,mppti) = Pwcs(t,mppti) + sum((1-worstcaseShFj).*Pjk,2);
                % Plsh(t,mppti) = Plsh(t,mppti) + sum(Pjk.*(1 - ShFj + kd.*ShFj),2);
                           
                if ~tricky(str)
                % Strings with uniform irradiance
                
                    strIVpp(str) = scale(ModIV.getIVpp(Gs,Ts),modulesperstring(mppti),1);

                    if opt.bypassdiodes || opt.modulediodes
                    % Add resulting IV-curve of all diodes in a string
                        strbpIVpp = scale(Diode.getIVpp(Ts),diodespermodule(1)*modulesperstring(mppti),diodespermodule(2),arrIVLimits);
                        strIVpp(str) = addparallel([strbpIVpp,strIVpp(str)]);
                    end
                    if opt.blockingdiodes
                    % Add a diode in series, if required
                        strIVpp(str) = addseries([Diode.getIVpp(Tarr_avg),strIVpp(str)]);
                    end
                    simtimer.add2lap('string');
                    continue;
                end
                
                % indices of partly-shaded modules, only shaded if there is a beam component!
                ispartshaded = ShFj > 0 & ShFj < 1;  
                ispartshaded = ispartshaded & ((Gush - Gsh)./Gavg > opt.Gvartol);
                if opt.celldefects
                    ispartshaded = ispartshaded | defective(in_str)';
                end

                % Classify modules without beam shading in clusters of approx. equal irradiance
                % the idea is to reduce number of unique IV curves and thus the computational expense
                % Clustering (dg/g)² < RelTol ensures that Ns·Pmp(Gavg)/sum(Pmp(G)) - 1 ~ RelTol
                curvetypes = zeros(nm,1);
                if nnz(~ispartshaded) > 1
                    curvetypes(~ispartshaded) = clusterdata(Gavg(~ispartshaded)'/mean(Gavg),...
                        'cutoff',opt.Gvartol,'criterion','distance','linkage','complete');
                else
                   curvetypes(~ispartshaded) = 1; 
                end

                % Use negative indices for beam-shaded modules (force calculation of individual IV curve)
                curvetypes(ispartshaded) = -(1:nnz(ispartshaded));

                simtimer.add2lap('string');

                modIVpp(:) = mdpwl();

                % Generate IV curves for modules under uniform-irradiance (positive stamps)
                npieces = 0;
                for s = 1:max(curvetypes)
                    inSbin = (curvetypes == s);
                    Gs = mean(Gavg(inSbin));   % Average irradiance for all modules in the bin
                    Ts = mean(Tavg(inSbin));   % Average cell-temperature ...
                    pp = ModIV.getIVpp(Gs,Ts);

                    if opt.bypassdiodes || opt.modulediodes
                    % Add resulting IV-curve of all diodes in a module
                        modbpIVpp = scale(Diode.getIVpp(Ts),diodespermodule(1),diodespermodule(2),modIVLimits);
                        modIVpp(inSbin) = addparallel([pp,modbpIVpp],modIVLimits,0); % skip simplification (faster)
                    else
                      modIVpp(inSbin) = scale(pp,1,1,modIVLimits,0);  % just clip it to bounds
                    end
                    npieces = npieces + modIVpp(find(inSbin,1)).n;
                end

                % Calculate IV curves for defective/partially-shaded modules
                for s = find(curvetypes < 0)'
                    [tr,p] = ind2sub([Nu,Np],in_str(s));   % physical location of the module

                    Gs = double(POA.Dpoa(t,tr,p)+[0,POA.Bpoa(t,tr,p)]); % G(shaded,not-shaded)
                    Ts = double([POA.Tsh(t,tr,p),POA.Tush(t,tr,p)]);    % Tc(shaded,not-shaded)

                    if ispartshaded(s)
                    % Partially-shaded module

                        % Get shading polygon(s) for module p, and their weights (if interpolated)
                        pidx = find(ShRes.BLPidx(ShRes.partidx(t,tr),:));
                        blpoly = ShRes.BLPoly(pidx);
                        blpwt = full(ShRes.BLPidx(ShRes.partidx(t,tr),pidx));

                        %shf = ShadingFactors(MountGeom.elements(p),blpoly,opt.RelTol);
                        %shf = QuickShadingFactors(Xcc(p,:,:),Ycc(p,:,:),blpoly);
                        shf = CellShadingFactors(p,blpoly,blpwt,opt.quickshading);
                    else
                    % defective
                        shf = zeros(Nepm,Ncpe);
                    end

                    if defective(tr,p)
                        def = full(reshape(POA.Defects(Nu*(p-1)+tr,:),Nepm,Ncpe));
                        kd = Gs(1)/Gs(2);
                        % Calculate an equivalent shading-factor, as if Dpoa = 0
                        Gs(1) = 0;
                        shf = (1-kd)*(1-def).*shf + def;
                    end

                    % get custom PV curve when there's no other choice 
                    pp = shadedIVcurve(shf,Gs,Ts);
                    if opt.modulediodes, pp = addparallel(Diode,pp,Ts(2)); end
                    modIVpp(s) = pp;
                    npieces = npieces + pp.n;
                end
                simtimer.add2lap('curves');

                % PROVISIONAL: the PWL-simplification algorithm is currently so slow, that when
                % merging PWL-curves with less than O(10/RelTol) segments makes little sense. 
                % If it looks like this will be the case, skip simplification alltogether.
                % This might change when mdpwl.douglaspeucker is implemented as a MEX function
                if Ns*npieces < 10/opt.RelTol, mergetol = 0; else, mergetol = Tol.*[Nm 1 1]; end

                % Add all module-IV-curves in series
                strIVpp(str) = addseries(modIVpp,arrIVLimits,mergetol);
                
                if opt.blockingdiodes
                    % Add a diode in series, if required
                    strIVpp(str) = addseries([Diode.getIVpp(Tarr_avg),strIVpp(str)]);
                end
                simtimer.add2lap('string');
                
                if Psh(t,mppti,str) > Push(t,mppti,str)
                    err = (Psh(t,mppti,str)/Push(t,mppti,str)-1)*100;
                    if abs(err) > opt.RelTol
                        warning('TimeSeriesCalc:PshPush_str',...
                            'Psh > Push by %0.2f%% at t:%d, i:%d, s:%d',err,t,mppti,str); 
                    else
                        Push(t,mppti,str) = Psh(t,mppti,str);
                    end
                    if ~isempty(dbstatus()), keyboard(); end % ### DEBUG!
                end
            end

            % when done with all strings, add them in parallel - do not simplify resulting curve, yet (*)
            arrIVpp = addparallel(strIVpp(used),arrIVLimits,0); %...,mergetol);

            % Find Vmp and check for inverter MPPT limitations
            [~,Vdc(t,mppti),~,Vsh(t,mppti)] = arrIVpp.mpp(Vlim);

            simtimer.add2lap('array');

            if opt.savingIVpps || plotobj.plotting
            % % remember a simplified (within Tol) version of the IV curve (*)
                % arrIVpp = mdpwl(arrIVpp.x,arrIVpp.y,Tol.*[Nm Ns 1]);
                if opt.savingIVpps, savedIVpps(t,mppti) = arrIVpp; end
                simtimer.add2lap('plot');
            end

            % Go back to the string curves and get individual ouputs
            for str = find(used)
                Ik = strIVpp(str).val([Vdc(t,mppti),Vsh(t,mppti)]);
                Pdc(t,mppti,str) = Ik(1)*Vdc(t,mppti);
                Psh(t,mppti,str) = Ik(2)*Vsh(t,mppti);
            end
            simtimer.add2lap('string');

            if sum(Psh(t,mppti,:)) > sum(Push(t,mppti,:))
                err = (sum(Psh(t,mppti,:))/sum(Push(t,mppti,:))-1)*100;
                if abs(err) > opt.RelTol
                    warning('TimeSeriesCalc:PshPush','Psh > Push by %0.2f%% at t:%d, i:%d',err,t,mppti); 
                end
                if ~isempty(dbstatus()), keyboard(); end % ### DEBUG! 
            end
            simtimer.add2lap('array');

            % every once in a while keep track of the time
            progresscount(1) = progresscount(1) + nnz(tricky);
            if ~plotobj.plotting && plotobj.worthupdating
                
                simprogress = (mppti+(t-t0)*Ni)/((Nt-t0+1)*Ni); % all-timesteps
                simprogress = 0.2*simprogress + 0.8*progresscount(1)/progresscount(2);
                progressmsg = sprintf('Inverter %d/%d, Time-Step %d/%d',mppti,Ni,t,Nt);
                plotobj.updatewaitbar(simprogress,progressmsg,'-addtime');
                
            elseif plotobj.plotting
                plotobj.reset();
                
                colors = autumn(numel(strIVpp));
                colors(:,4) = 0.2;
                for j = find(used)
                    plotobj.plotarraycurves(scale(strIVpp(j),1,stringspermppt(mppti)),'color',colors(j,:));
                    % [p,v] = strIVpp(j).mpp();
                    % plotter.plotpointoncurves(v,p*stringspermppt(mppti),'x','color',colors(j,:));
                    % plotter.plotpointoncurves(Vdc(t,mppti),,'ro');
                end
                
                [Gbc,idx] = max(double(POA.Dpoa(t,allmodules)+POA.Bpoa(t,allmodules)));
                Tbc = POA.Tush(t,allmodules(idx));
                pp = scale(ModIV.getIVpp(Gbc,Tbc),modulesperstring(mppti),stringspermppt(mppti));
                plotobj.plotarraycurves(pp,'g:');
                pp = scale(ModIV.getIVpp(Gwc,Twc),modulesperstring(mppti),stringspermppt(mppti));
                plotobj.plotarraycurves(pp,'r:');
                plotobj.plotarraycurves(arrIVpp,'b-');
                plotobj.plotpointoncurves(mean(Vush(t,mppti,:),'omitnan'),sum(Push(t,mppti,:),'omitnan'),'go');
                plotobj.plotpointoncurves(Vdc(t,mppti),sum(Pdc(t,mppti,:),'omitnan'),'mo');
                plotobj.plotpointoncurves(Vsh(t,mppti),sum(Psh(t,mppti,:),'omitnan'),'bo');

                % tridx = unique(ArrDef.psubs(mppti,:,:)*[1;0]); % index of trackers in mppti
                % shwt = cell(1,numel(tridx));
                % shp = cell(1,numel(tridx));
                % for j = 1:numel(tridx)
                %     if ShRes.fullshaded(t,tridx(j))
                %         shp{j} = -1;
                %         shwt{j} = 1;
                %     elseif ShRes.partshaded(t,tridx(j))
                %         pidx = find(ShRes.BLPidx(ShRes.partidx(t,tridx(j)),:));
                %         shwt{j} = full(ShRes.BLPidx(ShRes.partidx(t,tridx(j)),pidx));
                %         shp{j} = ShRes.BLPoly(pidx);
                %     else
                %         shp{j} = 1;
                %         shwt{j} = 1;
                %     end
                % end
                % plotter.plotshadedtrackers(shp,shwt);

                plotobj.addtextline(sprintf('t: %d/%d (%s)',t,Nt,timelabels(t,:)));
                plotobj.addtextline(stepinfo{t});
                plotobj.addtextline(sprintf('mppt-i: %d/%d',mppti,Ni));
                plotobj.addtextline('');
                plotobj.addtextline(sprintf('Push: %0.2f kW',sum(Push(t,mppti,:),'omitnan')/1000));
                plotobj.addtextline(sprintf('Plsh: %0.2f kW',Plsh(t,mppti)/1000));
                plotobj.addtextline(sprintf('Psh: %0.2f kW',sum(Psh(t,mppti,:),'omitnan')/1000));
                plotobj.addtextline(sprintf('Pdc: %0.2f kW',sum(Pdc(t,mppti,:),'omitnan')/1000));
                plotobj.addtextline(sprintf('Pwcs: %0.2f kW',Pwcs(t,mppti)/1000));
                plotobj.addtextline('');
                plotobj.addtextline(sprintf('Vush(max): %0.2f V',max(Vush(t,mppti,:))));
                plotobj.addtextline(sprintf('Vsh: %0.2f V',Vsh(t,mppti)));
                plotobj.addtextline(sprintf('Vdc: %0.2f V',Vdc(t,mppti)));
                plotobj.addtextline('');
                plotobj.addtextline(sprintf('strength: %0.1f%%',...
                    100*(Plsh(t,mppti)- sum(Psh(t,mppti,:),'omitnan'))/(Plsh(t,mppti)-Pwcs(t,mppti))));

                drawnow()
                if plotobj.pausing, plotobj.pause(); end
                if opt.exportplots
                    if size(timelabels,2)>1
                        figname = timelabels(t,[1:4 6:7 9:10 12:13 15:16]);
                    else
                        figname = sprintf('t%04d',t);
                    end
                    figname = sprintf('%s_mppt%02d.png',figname,mppti);
                    plotobj.exportfigure(figname);
                end
            end
            plotobj.refresh(); % apply any queued request to change plotter instance        
        end
        simtimer.add2lap('plot');

        if opt.solverlog
        % Step-by-step output Log

          if ~ismember(pvsLog,fopen('all')), pvsLog = fopen('pvSolverLog.txt','a'); end 

          fprintf(pvsLog,'\r\nTimestep %d/%d (%s,%s) - (MPPTi): Value',t,Nt,timelabels(t,:),stepinfo{t});
          fprintf(pvsLog,'\r\n\tPush [W]:\t');
          fprintf(pvsLog,'(%02d):%0.2f\t',[1:Ni;sum(Push(t,:,:),3)]);
          fprintf(pvsLog,'\r\n\tPwcs [W]:\t');
          fprintf(pvsLog,'(%02d):%0.2f\t',[1:Ni;Pwcs(t,:)]);
          fprintf(pvsLog,'\r\n\tPsh  [W]:\t');
          fprintf(pvsLog,'(%02d):%0.2f\t',[1:Ni;sum(Psh(t,:,:),3)]);
          fprintf(pvsLog,'\r\n\tPlsh [W]:\t');
          fprintf(pvsLog,'(%02d):%0.2f\t',[1:Ni;Plsh(t,:)]);
          fprintf(pvsLog,'\r\n\tPdc  [W]:\t');
          fprintf(pvsLog,'(%02d):%0.2f\t',[1:Ni;sum(Pdc(t,:,:),3)]);
          fprintf(pvsLog,'\r\n\tVush(avg) [V]:\t');
          fprintf(pvsLog,'(%02d):%0.2f\t',[1:Ni;mean(Vush(t,:,:),3)]);
          fprintf(pvsLog,'\r\n\tVsh  [V]:\t');
          fprintf(pvsLog,'(%02d):%0.2f\t',[1:Ni;Vsh(t,:,:)]);
          fprintf(pvsLog,'\r\n\tVdc [V]:\t');
          fprintf(pvsLog,'(%02d):%0.2f\t',[1:Ni;Vdc(t,:,:)]); 
          fprintf(pvsLog,'\r\n');

          simtimer.add2lap('log');
        end

        % Dump complete workspace to partial-result file every AUTOSAVEMIN minutes
        if ~isempty(opt.backup) && simtimer.globaltime - lastsaved > opt.autosavemin*60 || t == Nt
            save(opt.backup,VARS_TO_SAVE{:});       
            lastsaved = simtimer.globaltime;
            simtimer.add2lap('backup');
        end

    end
    delete(cleanupobj); % plotter.close(true);

    SolRes.Pdc = Pdc;
    SolRes.Vdc = Vdc;
    SolRes.Psh = Psh;
    SolRes.Vsh = Vsh;
    SolRes.Push = Push;
    SolRes.Vush = Vush;
    SolRes.Plsh = Plsh;
    SolRes.Pwcs = Pwcs;
    SolRes.simtimer = simtimer;
    if opt.savingIVpps, SolRes.arrIVpp = savedIVpps; end

    if opt.solverlog, fclose(pvsLog); end
    simtimer.add2lap('packing');

    % end of main function

    function IVpp = shadedIVcurve(shf,gc,tc)
    % Return the IV curve of a PV-module whose cells are shaded according to shf
    % cellint - cell IV curve ModuleInterpolant object
    % shf - Ns·Nc array of shading factors (0 = not shaded, 1 = shaded)
    % gc - POA irradiance, gc(1) = diffuse, gc(2) = global
    % tc - Cell temperature, tc(1) = shaded cell, tc(2) = unshaded cell
    % opt - structure with fields:
    %   .cellshadingtol (smallest difference to consider between two shaded cells)
    %   .bypassdiodes (bypass diodes in every cell-block)
    %   .parallelelements (cell-blocks connected in parallel, not series)
    %   .modulediodes (bypass diode for complete module)
    % blckBPDpp - PWL (pp-spline) object, required if opt.bypassdiodes or modulediodes

        if opt.bypassdiodes || opt.modulediodes
        % Get bypass-diode IV-curve assuming package temp. ~ unshaded-cell temp.
            blckBPDpp = getIVpp(Diode,tc(2));
            blckBPDpp = scale(blckBPDpp,1,1,modIVLimits); 
        end
        % if there are no bypass diodes, and substrings are in series, then
        % we can treat the module as a single cell-substring
        if ~opt.bypassdiodes && ~opt.parallelelements, shf = shf(:); end

        gc = gc(2)-(gc(2)-gc(1))*shf';
        tc = tc(2)+(tc(2)-tc(1))*shf';

        % Classify cells in clusters of approx. equal irradiance
        % the idea is to reduce number of unique IV curves and thus the computational expense
        % Clustering (dg/g)² < RelTol ensures that N·Pmp(Gavg)/sum(Pmp(G)) - 1 ~ RelTol
        ctol = min(opt.cellshadingtol,opt.Gvartol);
        ctypes = clusterdata(gc(:)/mean(gc(:)),'cutoff',ctol,'criterion','distance','linkage','complete');
                
        if max(ctypes) == 1
        % Shortcut: No real mismatch in complete module: e.g. overcast sky
            IVpp = getIVpp(ModIV,min(gc(:)),mean(tc(:))); % use MIN irradiance!
            if opt.bypassdiodes || opt.modulediodes
                % Add resulting IV-curve of all diodes in a module:
                blckBPDpp = scale(blckBPDpp,diodespermodule(1),diodespermodule(2),modIVLimits);
                IVpp = addparallel([IVpp,blckBPDpp],modIVLimits,0); % skip simplification
            else
                IVpp = scale(IVpp,1,1,modIVLimits,0);  % just clip it to bounds
            end
            return;
        end

        % Evaluate the IVpps for each cell-state
        ivpps(max(ctypes),1) = mdpwl();
        for u = 1:max(ctypes)
            inu = ctypes == u;
            ivpps(u) = getIVpp(cellIV,min(gc(inu)),mean(tc(inu))); % use MIN irradiance!
        end

        [ns,nc] = size(shf); % cell-substrings-per-module, cells-per-substring
        ctypes = reshape(ctypes,nc,ns);

        % add-up the cell-IVpps for every substring
        subsIV(ns,1) = mdpwl();
        for u = 1:ns
            if all(ctypes(:,u)==ctypes(1,u))
            % No mismatch on sub-string
                subsIV(u) = scale(ivpps(ctypes(1,u)),nc,1,modIVLimits,0);
            else
            % Cell-wise series connection
                subsIV(u) = addseries(ivpps(ctypes(:,u)),modIVLimits,0);
            end
            % add bypass diode
            if opt.bypassdiodes, subsIV(u) = addparallel([subsIV(u),blckBPDpp],modIVLimits,0); end
        end

        % add-up substrings in module
        if opt.parallelelements, IVpp = addparallel(subsIV,modIVLimits,0);
        else, IVpp = addseries(subsIV,modIVLimits,0);
        end
        if opt.modulediodes, IVpp = addparallel([IVpp,blckBPDpp],modIVLimits,0); end
    end
end

function [F,Nepm,Ncpe] = CellShadingFactors(p,blp,wt,quickshading)
% [~,NEPM,NCPE] = CELLSHADINGFACTORS(PCKGEOM,Np) - Parse packed-pvArea PCKGEOM, and initialize 
%   persistent copies of cell-vertices ([Ne,Nc,m,Np] arrays), for efficient call to INPACKEDPOLYGON 
%   during runtime.
%
% F = CELLSHADINGFACTORS(P,BLP,W,PCKGEOM,QUICK) - Return an [Nepm, Ncpe] array of factors from 0 
%   (no shade) to 1 (fully shaded cell), for the Nepm·Ncpe cells of module P in geometry PCKGEOM,
%   [pack16'd pvArea] applying light polygon(s) BLP with weights WT. 
%
%   For step-by-step shading calculations, BLP = {Q} should contain a single polygon Q, and WT 
%   should be a scalar 1. For interpolated shading results, BLP = {A,B,..} should typically contain
%   1 to 3 polygons evaluated at nearby solar positions, with WP providing the weight of each
%   position (e.g. as the result of INTERPMATRIX(..,'-sph').

    persistent Xcc Ycc Pcc cellarea
    
    if nargin == 2 || isempty(cellarea)
    % A = CELLSHADINGFACTORS(PCKGEOM,Np)
    
        try
            assert(isstruct(p) && all(isfield(p,{'border','elements','dims','depth'})) && ...
                p.depth == 4 && p.dims(1) == blp && all(isfield(p.border,{'x','y'})) && ...
                isa(p.border.x,'int16'));

            Np = p.dims(1);   % Modules per Tracker
            Nepm = p.dims(2); % Elements-per-Module
            Ncpe = p.dims(3); % Cells-per-Element
            Ncv = numel(p.elements(1).elements(1).elements(1).border.x); % first cell's numel(x)

            Pcc = pack16(polygon());
            Pcc(Np,Nepm,Ncpe) = pack16(polygon());
            cellarea = zeros(Np,Nepm,Ncpe);
            Xcc = zeros(Nepm,Ncpe,Ncv,Np);
            Ycc = zeros(Nepm,Ncpe,Ncv,Np);
            for i = 1:Np
                for j = 1:Nepm
                    for k = 1:Ncpe
                        Pcc(i,j,k) = p.elements(i).elements(j).elements(k).border;
                        cellarea(i,j,k) = polygon.packedarea(Pcc(i,j,k));
                        Xcc(j,k,:,i) = double(Pcc(i,j,k).x)/double(Pcc(i,j,k).scale);
                        Ycc(j,k,:,i) = double(Pcc(i,j,k).y)/double(Pcc(i,j,k).scale);
                    end
                end
            end
            
            F = sum(cellarea,'all')/polygon.packedarea(p.border);
            assert(F <= 1,'ShRes.geom.border area is less than the cell of its cells!');
            if F < 0.75
                warning('Cell-area/Mount-area (%0.1f%%) seems exceptionally low',F*100);
            end

            F = []; 
            return;
        catch ERR
            throwAsCaller(MException('pvArraySolver:mountgeom',...
                ['Failed to initialize CellShadingFactors. Check that ShRes.mountgeom is ',...
                 'a packed, %d-element, depth-4 pvArea object: %s'],blp,getReport(ERR)));
        end
    else
        [Nepm,Ncpe,Ncv,~] = size(Xcc);
    end
    
    narginchk(4,4);
    
    % DEBUG
    if numel(blp)~=numel(wt) || abs(sum(wt)-1) > 1e-6 
    	if ~isempty(dbstatus()), keyboard(); end
    end

    F = zeros(Nepm,Ncpe);
    for k = 1:numel(blp)
        
        % Approx. factors, based on the number of cell-border vertices that are shaded
        f = 1-sum(polygon.inpackedpolygon(blp{k},Xcc(:,:,:,p),Ycc(:,:,:,p)),3)/Ncv;
        
        if ~quickshading
            partshadedcell = find(f > 0 & f < 1); % index of partially-shaded cells
            if any(partshadedcell)
            % Refine partially-shaded-cell factors
                [e,c] = ind2sub([Nepm,Ncpe],partshadedcell);
                for h = 1:numel(e)
                    a = polygon.packedarea(polyclip(Pcc(p,e(h),c(h)),blp{k},1));
                    f(partshadedcell(h)) = 1-a/cellarea(p,e(h),c(h));
                end
            end
        end
        F = F + f*wt(k);
    end
end

% function [X,Y] = centerpoints(obj)
% 	if isempty(obj.elements)
% 		[~,~,c] = rectangleproperties(obj.border);
% 		X = c(1); Y = c(2);
% 	else
% 		s = obj.dims;
% 		if numel(s) == 1, s = [s,1]; end
% 		X = zeros(s);
% 		Y = zeros(s);
% 		for j = 1:numel(obj.elements)
% 			[xj,yj] = centerpoints(obj.elements(j));
%             try
%                 X(j,:) = xj(:)';
%                 Y(j,:) = yj(:)';
%             catch
%                 keyboard;
%             end
% 		end
% 	end
% end

