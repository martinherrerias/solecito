function ok = checkdesign(varargin)
% OK = CHECKDESIGN(Trackers,ArrayDef,ModIVint,celltemp,Inverter,MeteoData,SunPos) - Perform a 
%   simplified, reduced (random sample) back-of-the-(virtual)envelope calculation to make sure 
%   that meteo-data, array-definition, and choice of module (MOD) and INVERTER all make sense.  
%   That is: from a system design's perspective -- individual arguments are expected to have been 
%   checked already for correctness and strict compatibility (CHECKARRAYDEFINITION, CHECKODM, ...)
%
%   More specifically, inverter losses (as returned by ONDEVAL) are expected to fall below 
%   predefined thresholds (5% default for most forms of clipping: clip_i, clip_v, clip_t, and
%   and thresh; 10% for clip_p; and 15% due to efficiency). Specific limits can be altered by
%   name,value pairs, e.g. 'clip_p',0.2.
%
%   All arguments can be omitted (or left empty), in which case they will be searched for in the 
%   base workspace (case-sensitive names given above).
%
%   If METEODATA or SUNPOS are not available, checks will be performed based on 'common' operating 
%   conditions (400/800 W/m² and 20°C), otherwise POAIRRADIANCE with IAM = 1 and no shading for
%   a subset of points and mounts will be used. (Defaults: 'MaxPts',5000,'MaxTrck',20)
%
% OUTPUT: returns a boolean scalar (true if all checks are passed). Additionally, with a 
%   '-verbose' flag, prints out 95% bounds for Voc, Isc, Vmp, Pmp, and inverter efficiency.
%   Use '-soft' to catch any error and return it as a warning, with the resulting OK = false.
%
% FUTURE: could be used to unify consistency checks, similar to how PRINTSUMMARY seeks to unify
%   display. Taking over consistency checks of e.g. ARRAYDEF vs TRACKERS
% 
% See also: ONDEVAL, CHECKARRAYDEFINITION, CHECKODM

    ARGS = {'Trackers','ArrayDef','ModIVint','celltemp','Inverter','MeteoData','SunPos'};
    REQ = [1,1,1,1,1,0,0];

    [opt,varargin] = getflagoptions(varargin,{'-verbose','-soft'});
    [opt,varargin] = getpairedoptions(varargin,opt);
    if opt.soft
        try
            ok = checkdesign(varargin{:},'verbose',opt.verbose);
        catch ERR
            warning('checkdesign:softcrash','Design-check crashed, review %s: \n\n%s',...
                shortliststr(ARGS,[],10),getReport(ERR));
            ok = false;
        end
        return; 
    end
    
    opt.MaxPts = 5000;
    opt.MaxTrck = 20;
    
    % Define field-names, tags, and threshold-criteria for ONDEVAL's LOSS structure
    FLD = {'thresh', 'min. Power threshold clipping';
           'self', 'night/Aux. losses';
           'clip_v', 'MPPT voltage clipping';
           'clip_p', 'max. Power clipping';
           'clip_i', 'input current clipping';
           'clip_t', 'temperature derate clipping';
           'eff', 'efficiency losses'};

    opt.thresh = 0.05;
    opt.self = 0.05;
    opt.clip_i = 0.05;
    opt.clip_v = 0.05;
    opt.clip_p = 0.1;
    opt.clip_t = 0.05;
    opt.eff = 0.15;
    
    [opt,varargin] = getpairedoptions(varargin,opt);
    varargin(end+1:7) = {[]};
    assert(numel(varargin) == 7,'Unrecognized arguments');
    [Trackers,ArrayDef,ModIVint,celltemp,Inverter,MeteoData,SunPos] = deal(varargin{:});
    
    % Load missing arguments from base workspace
    missing = false(1,numel(ARGS));
    for j = 1:numel(ARGS)
        if isempty(eval(ARGS{j})), missing(j) = ~getfrombase(ARGS{j}); end
    end
    missing(~REQ) = false;
    if any(missing)
        error('checkdesign:models',shortliststr(ARGS(missing),'Failed to load argument','colon',':'));
    end
    
    printif = @(varargin) opt.verbose && fprintf(varargin{:});
    printif('\nPerforming design-check...\n');
    ok = true;
    
    idx = ArrayDef.esubs(:);
    MPS = accumarray(idx(:,1:2),1); % Ni x max(Ns) array of modules-per-string
    Ns = sum(MPS > 0,2);            % number of strings-per-inverter 
    Nm = max(MPS,[],2);             % ... and modules-per-string
    assert(all(MPS == Nm | MPS == 0,'all'),'Strings with different number of modules!');
    
    % Keep only unique combinations of Ns, Nm
    [~,idx] = unique([Ns,Nm],'rows');
    Nm = Nm(idx); Ns = Ns(idx);
    
    Plim = [Inverter.Ps0 Inverter.Pdc0];
    Vlim = [Inverter.MPPTLow Inverter.MPPTHi];
        
    % Try to find MeteoData
    if ~isempty(MeteoData)
        [Gpoa,Dpoa,MeteoData] = sample_irradiance(MeteoData,SunPos,Trackers,opt);
    
        % Get cell temperatures for shaded and not shaded elements
        [Nt,Ntr,Np] = size(Gpoa);
        Tush = zeros(Nt,Ntr,Np);
        Tsh = zeros(Nt,Ntr,Np);
        ta = single(repmat(MeteoData.Ta,Ntr*Np,1)); 
        vw = single(repmat(MeteoData.vw,Ntr*Np,1));
            Tush(:) = celltemp(Gpoa(:),ta,vw);  
            Tsh(:) = celltemp(Dpoa(:),ta,vw);
            
        % Check interpolant limits
        Glim = [min(Gpoa(:)),max(Gpoa(:))];
        Tlim = [min([Tush(:);Tsh(:)]),max([Tush(:);Tsh(:)])];
        if Glim(1) < ModIVint.Gref(1) || Glim(2) > ModIVint.Gref(end) || ...
           Tlim(1) < ModIVint.Tref(1) || Tlim(2) > ModIVint.Tref(end)
            warning('checkdesign:modIVlims',['Irradiance and/or cell-temperature values ',...
                '(%0.1f - %0.1f W/m², %0.1f - %0.1f °C) exceed IV-curve interpolant limits ',...
                '(%0.1f - %0.1f W/m², %0.1f - %0.1f °C). Check your input meteo-data & ',...
                'cell-temperature model, or rebuild your IV-curve interpolant with appropriate ',...
                'SimOptions.Gegrid / Tcgrid settings'],Glim(1),Glim(2),Tlim(1),Tlim(2),...
                ModIVint.Gref(1),ModIVint.Gref(end),ModIVint.Tref(1),ModIVint.Tref(end));
            ok = false;
        end
        
        % Check if inverter selection makes any sense
        notdark = Gpoa > getSimOption('minGHI');
        [Pmp,Vmp] = getMPP(ModIVint,Gpoa(notdark),Tush(notdark));
        
        % Get approx. Voc on shaded (colder, higher) conditions
        Voc = ModIVint.getVoc(Gpoa(notdark),Tsh(notdark)).*Nm';
        Voc = max(Voc,[],2);
        [~,idx] = max(Voc(:));
        printif('\tMax. string Voc is %0.1f V at %0.0f W/m² and %0.1f °C\n',Voc(idx),Gpoa(idx),Tsh(idx));
        
        % Get approx. Isc under no-shade conditions
        Isc = ModIVint.getIsc(Gpoa(notdark),Tush(notdark)).*Ns';
        Isc = max(Isc,[],2);
        [~,idx] = max(Isc(:));
        printif('\tMax. array Isc is %0.1f A at %0.0f W/m² and %0.1f °C\n',Isc(idx),Gpoa(idx),Tsh(idx));
  
        % Check Pmp and Vmp...
        Pmp = Pmp.*(Ns.*Nm)';
        Vmp = Vmp.*Nm';
        [Pac,Loss] = ONDeval(Inverter,Pmp,Vmp,ta(notdark));
        eta = Pac./Pmp;
        
        printif('\t95%% bounds for Pmp(not shaded) are %0.2fkW to %0.2fkW, inverter has %0.2fkW\n',...
                    prctilew(Pmp,2.5,Pmp)/1000, prctilew(Pmp,97.5,Pmp)/1000,Inverter.Pdc0/1000);
        printif('\t95%% bounds for Vmp(not shaded) are %0.1fV to %0.1fV, inverter has %0.1fV to %0.1fV\n',...
                    prctilew(Vmp,2.5,Pmp), prctilew(Vmp,97.5,Pmp),Inverter.MPPTLow,Inverter.MPPTHi);
        printif('\t95%% bounds for final inverter efficiency %0.1f%% to %0.1f%%\n',...
                    prctilew(eta,2.5,Pmp)*100, prctilew(eta,97.5,Pmp)*100);
        
        for j = 1:size(FLD,1)
            f = sum(Loss.(FLD{j,1}))/sum(Pac);
            if f > opt.(FLD{j,1})
                warning('checkdesign:losses','Inverter %s of %0.1f%% ',FLD{j,2},f*100);
                ok = false;
            end
        end
        
        if ok
            printif('\tElectrical models seem appropriate for meteo-data and array-definition.\n');
        end
    else
        % Perform common-sense checks only
        loadfactor = (ModIVint.Vmp0*ModIVint.Imp0*Ns.*Nm)/Inverter.PacMax;
 
        Tc = celltemp([400,800],[20,20],[1,1]);
        [Pmp,Vmp] = getMPP(ModIVint,[400,800],Tc);
        Pmp = Pmp.*Ns.*Nm;
        Vmp = Vmp.*Nm;
        eta = ONDeval(Inverter,Pmp,Vmp,20)./Pmp;
        
        if ~isscalar(Ns)
            Pmp = minmax(Pmp(:)');
            Vmp = minmax(Vmp(:)');
            eta = minmax(eta(:)');
            loadfactor = minmax(loadfactor(:)');
        end
        
        printif(['@ 400/800 W/m² and 20°C:\n',...
                '\tVmp is %0.1f V to %0.1f V, (inverter takes %0.1f V to %0.1f V)\n',...
                '\tPmp is %0.1f kW to %0.1f kW, (inverter is rated %0.1f kW)\n',...
                '\tInverter efficiency is %0.1f%% / %0.1f%%\n'],...
                Vmp(1),Vmp(2),Inverter.MPPTLow,Inverter.MPPTHi,...
                Pmp(1)/1000,Pmp(2)/1000,Inverter.Pdc0/1000,eta(1)*100,eta(2)*100);
        
        msg = {};
        if any(Pmp./(Ns.*Nm)./(ModIVint.area*[400,800]) < 0.1)
            msg{end+1} = 'module efficiency seems exceptionally low';
        end
        if any(Vmp < Vlim(1) | Vmp > Vlim(2) | Pmp < Plim(1) | Pmp > Plim(2))
            msg{end+1} = 'array''s MPP is outside inverter''s operating window';
        end
        if any(loadfactor < 0.9 | loadfactor > 1.5)
           msg{end+1} = sprintf('inverter load factor (%0.2f-%0.2f) is not normal',loadfactor);
        end
        if any(eta < 1-opt.eff)
            msg{end+1} = 'inverter efficiency seems abnormally low';
        end
        if ~isempty(msg)
            warning('checkdesign:pmpatnoc',['Check your system definition: ' shortliststr(msg)]);
            ok = false;
        end
        
        if ok, printif('\tDesign seems ok at standard-conditions\n'); end
    end
end

function gotitallright = getfrombase(varname,localvarname)
% Try to find the variable 'varname' in the base workspace, and assign it to 'localvarname' in the
% caller workspace.

    if nargin < 2, localvarname = varname; end
    
    if evalin('base',sprintf('exist(''%s'',''var'');',varname))
        assignin('caller',localvarname,evalin('base',varname));
        gotitallright = true;
    else
        %assignin('caller',varname,[]);
        gotitallright = false;
    end
end

function [Gpoa,Dpoa,MeteoData] = sample_irradiance(MeteoData,SunPos,Trackers,opt)
% Return (possibly a reduced sample of) unshaded plane of array irradiance (global Gpoa and 
% diffuse Dpoa), or, if that's not possible, simply Gpoa = GHI and Dpoa = DHI. 

    printif = @(varargin) opt.verbose && fprintf(varargin{:});

    F = MeteoData.GHI > 0;
    Nt = nnz(F);
    if Nt > opt.MaxPts
        F = find(F);
        F = F(randperm(Nt,opt.MaxPts));
    end
    Nt = numel(MeteoData.GHI);
    MeteoData = filterstructure(MeteoData,F,Nt);
    
    Gpoa = [];
    if ~isempty(SunPos)
        
        try
            SunPos = filterstructure(SunPos,F,Nt);

            Ntr = numel(Trackers.analysedtrackers);
            if Ntr > opt.MaxTrck
                F = randperm(Ntr,opt.MaxTrck);
                Trackers.analysedtrackers = Trackers.analysedtrackers(F);
            end

            POA = poairradiance(MeteoData,SunPos,Trackers,[],[],'fIAM',1);
            Dpoa = single(POA.Dpoa0); 
            Gpoa = single(POA.Bpoa0) + Dpoa;
            printif('Using unshaded POA-irradiance\n');
        catch ERR
            warning('Failed to calculate unshaded POA irradiance:\n\n%s', getReport(ERR));
        end
    end

   if isempty(Gpoa)  % the best we can do
       printif('Using horizontal irradiance\n');
       Gpoa = single(MeteoData.GHI); 
       Dpoa = single(MeteoData.DHI);
   end 
end
