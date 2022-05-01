function [ctf,ctf0,C] = getcelltempfcn(varargin)
% [CTF,ctf0,C] = GETCELLTEMPFCN(S) - Retrieve cell-temperature model parameters from structure
%   S. Check values, complete with defaults, and create cell-temperature-function handles
%   CTF, CTF0 for MPP and VOC cell temperatures, respectively -- both of the form @(Ge,Ta,Vw).
%   C is the complete set of parameters.
%
% [CTF,ctf0,C] = GETCELLTEMPFCN(MDL,['par',val]) - use defaults for a given model, modified
%   by name-value pairs. Equivalent to S = struct('model',MDL,'par',val,..).
%
% [CTF,ctf0,C] = GETCELLTEMPFCN() - load defaults from getSimOption('CellTemp').
%
% Supported models (and their default parameters) are given below.
% NOTE: setSimOption('CellTemp',struct('par',val,..)) can be used to override these defaults.
% Priority is: name-value arguments > first structure argument > SimOptions > defaults.
%
% Sandia Cell Temp. model [1] (see PVL_SAPMCELLTEMP):
%
%     S.model = 'SAPM'
%     S.type = 'Glass/cell/glass - Open rack'*
%     S.a_wind = -3.47
%     S.b_wind = -.0594	
%     S.delT = 3.0 [K]
%
%     (*) GETCELLTEMPFCN('sapm',[type]) - [Interactive] pick of SAPM temperature model defaults,
%        based on module materials and mounting types.
%
% Faiman quasi-physical model [2] (see CELLTEMPUVALUES):
%
%     S.model = 'Uval'
%     S.Uconst = 29.0   Heat transfer coefficient (wind independent) [W/m²K]
%     S.Uwind = 0.0     Heat transfer coefficient (wind dependent) [W·s/m³K]
%     S.absort = 0.9    Optical absorptivity
%     S.eff = 0.16 (§)  Module efficiency (function of irradiance and cell-temp.)
% 
%     (§) GETCELLTEMPFCN(ODM,'model','uval',...) - where ODM is a valid One-Diode-Model
%       parameter-set (see CHECKODM), uses a default efficiency function:
%
%       S.eff = @(G,Tc) onediodeMPP(translateODM(ODM,G,Tc))./(ODM.area*G); 
%
%     Faiman model with Uconst adjusted to NOCT (Uconst has precedence! -- i.e. must be unset)
%
%     S.NOCT = 44.8;     % Cell Temperature at NOC* (°C)
%     S.isCPV = false;   % (*) Adjusts NOC standard
% 
% Mattei et al. model [3] (see CELLTEMP_MATTEI):
%
%     S.model = 'mattei';
%     S.Uconst = 24.1;     % Intercept of heat transfer coefficient [W/m²K]
%     S.Uwind = 2.9;       % Slope of heat transfer coefficient [(W/m²K)/(m/s)]
%     S.absort = 0.81;	   % tau·alpha - Absortivity x IR Emissivity of module    
%     S.eta_ref = 0.16;    % efficiency at STC
%     S.beta = 0.004;      % Pmp temperature derate factor [1/K] (-muPmp in ODM)
%     S.gamma = 0.0;       % Pmp irradiance log. derate factor (typ. 0 - 0.12)
%
%     Mattei model with Uconst adjusted to NOCT (Uconst has precedence! -- i.e. must be unset)
%
%     S.NOCT = 46.4;     % Cell Temperature at NOC* (°C)
%     S.isCPV = false;   % (*) Adjusts NOC standard
%
% OUTPUT: 
%   CTF - cell temperature function @(g,t,v) of irradiance, ambient temperature, and wind speed.
%   CTF0 - cell temperature function @(g,t,v) at zero module efficiency (e.g. for NOCT)
%   C - [completed] parameter structure.
%
% References
%   [1] King, D. et al, 2004, "Sandia Photovoltaic Array Performance Model", SAND Report
%       3535, Sandia National Laboratories, Albuquerque, NM
%   [2] D. Faiman, "Assessing the outdoor operating temperature of photovoltaic modules,"
%       Progress in Photovoltaics: Research and Applications, vol. 16, no. 4, pp. 307–315, 2008.
%   [3] Mattei, M., Notton, G., Cristofari, C., Muselli, M., Poggi, P., 2006. Calculation   
%       of the polycrystalline PV module temperature using a simple method of energy balance. 
%       Renewable Energy 31, 553–567. https://doi.org/10.1016/j.renene.2005.03.010
%
% TODO: this should probably be a class?
%
% EXAMPLES:
%
% [CTF,ctf0,C] = getcelltempfcn('sapm','type','Glass/cell/glass - Open rack')
% [CTF,ctf0,C] = getcelltempfcn('noct',45,'uwind',2)
% ODM = checkODM(); [CTF,ctf0,C] = getcelltempfcn(ODM{1},'model','Uval')
%
% See also: PVL_SAPMCELLTEMP, CELLTEMPUVALUES, CELLTEMP_MATTEI

    ALIAS = {'Uval',{'uval','faiman','celltempuvalues'};
             'Mattei',{'mattei','celltemp_mattei'}
             'SAPM',{'sapm','pvl_sapmcelltemp','sandia','king'}};

    % Load SAPM Cell-Temperature Model parameters, as default (just in case)
    [DEF.SAPM(1:6).type] = deal('Glass/cell/glass - Open rack',...
                             'Glass/cell/glass - Close roof mount',...
                             'Glass/cell/polymer sheet - Open rack',...
                             'Glass/cell/polymer sheet - Insulated back',...
                             'Polymer/thin-film/steel - Open rack',...
                             '22X Linear Concentrator - Tracker');
    [DEF.SAPM(:).a_wind] = deal(-3.47,-2.98,-3.56,-2.81,-3.58,-3.23);
    [DEF.SAPM(:).b_wind] = deal(-.0594,-.0471,-.0750,-.0455,-.113,-.130);
    [DEF.SAPM(:).delT] = deal(3,1,3,0,3,13);
    [DEF.SAPM(:).isCPV] = deal(0,0,0,0,0,1);
    % NOCT = arrayfun(@(E,ws,ta,a,b,dT) round(pvl_sapmcelltemp(E,1000,a,b,ws,ta,dT),1),...
    %                 [800 800 800 800 800 900],[1 1 1 1 1 2],[20 20 20 20 20 20],...
    %                 [DEF.SAPM.a_wind],[DEF.SAPM.b_wind],[DEF.SAPM.delT],'unif',0);
    % [DEF.SAPM(:).NOCT] = deal(NOCT{:});
    
    DEF.Uval.type = '';
    DEF.Uval.Uconst = 29.0;
    DEF.Uval.Uwind = 0.0;
    DEF.Uval.absort = 0.9;
    DEF.Uval.eff = 0.16;
    DEF.Uval.isCPV = false;
    DEF.Uval.NOCT = celltempUvalues(800,20,1,DEF.Uval,0);

    [~,~,DEF.Mattei] = celltemp_mattei(0,0,0);
    DEF.Mattei.type = '';
    DEF.Mattei.isCPV = false;
    DEF.Mattei.NOCT = celltemp_mattei(880,20,1,'eta_ref',0);

    RELTOL = getSimOption('RelTol');
    
    % For (dP/dT)/P ~ -0.004 1/K and dT = 100·tol -> |dP/P| = dT·(dP/dT)/P ~ |0.4·tol|
    TOL_TEMP = RELTOL*100;
    
    KFIELDS = cellfun(@fieldnames,struct2cell(DEF),'unif',0);
    KFIELDS = unique(cat(1,{'model'}',KFIELDS{:}),'stable');
    
    opt =  getSimOption('CellTemp');
    if ischar(opt), opt = checkmodelname(struct('model',opt)); 
    else
        parsestruct(opt,{},'opt',setdiff(KFIELDS,{'model','type'}),'numeric','real','scalar');
        parsestruct(opt,{},'opt',{'model','type'},'class','char');
        opt = checkmodelname(opt,'-soft');
    end

    isODM = false;
    if isstruct(varargin{1})
        S0 = checkmodelname(varargin{1},'-soft');
        
        wr = naptime('all','off'); %#ok<NASGU>
        try %#ok<TRYNC>
            S0 = checkODM(S0);
            effODM = @(g,t) onediodeMPP(translateODM(S0,g,t),TOL_TEMP/4)./(S0.area*g);
            isODM = true; 
        end 
        clear wr
        
        varargin(1) = [];
    else
        S0 = struct(); 
    end
    nvp = getpairedoptions(varargin,KFIELDS,'dealrest',1); % name-value-pairs
    nvp = checkmodelname(nvp,'-soft');
    
    if ~isfield(nvp,'model')
        if isfield(nvp,'NOCT'),nvp.model = 'Mattei';
        elseif isfield(S0,'model')
            nvp.model = S0.model;
        else
            assert(isfield(opt,'model'),'Undefined cell-temperature model');
            nvp.model = opt.model;
        end
    end
    
    % remove irrelevant fields (with a warning)
    nvp = parsemdl(nvp,nvp.model,'provided arguments');
    opt = parsemdl(opt,nvp.model,'SimOptions.CellTemp');
    S = parsemdl(S0,nvp.model,'provided structure');
    
    % extract efficiency parameters from ODM
    if strcmp(S.model,'Mattei') && isODM
        [ok,P] = testefficiency(effODM);
        if ok, S = completestruct(S,P); end
    end
    
    % ... and merge with priority: nvp > S > opt
    C = completestruct(nvp,S);
    C = completestruct(C,opt);
    S = completestruct(S,opt);
    opt = completestruct(opt,S);
    [~,~,txt] = comparestruct(C,mergestructures(S,opt),'and');
    if ~isempty(txt)
        warning('Inconsistent field(s) in cell-temperature model:\n%s\n',txt); 
    end

    % Get/check efficiency function (for Faiman model)
    if strcmp(C.model,'Uval')
        if isfield(C,'eff')
            if isnumeric(C.eff)
                validateattributes(C.eff,'numeric',{'scalar','real','>=',0','<=',1});
            else
                assert(testefficiency(C.eff),'Invalid efficiency function');
            end
        elseif isODM
            try  %#ok<TRYNC>
                assert(testefficiency(effODM));
                C.eff = effODM;
            end
        end
    end

    % (interactive) SAPM model default pick
    if strcmpi(C.model,'SAPM')
        if all(isfield(C,setdiff(fieldnames(DEF.SAPM),{'type','isCPV'})))
            if ~isfield(C,'type'), C.type = 'Custom'; end
            n = 1;
        else
            if isfield(C,'type')
               [~,n] = parselist(C.type,{DEF.SAPM.type},'SAPM model type','');
            elseif runningfromUI()
                n = listdlg('PromptString','Choose module type and mounting','ListSize',[300,100],...
                        'ListString',{DEF.SAPM(:).type},...
                        'SelectionMode','single','Name','SAPM Cell-Temp. Model');
                assert(~isempty(n),'Invalid Module/mointing type!');
            else, n = 1;
            end
            C.type = DEF.SAPM(n).type;
        end
        DEF.SAPM = DEF.SAPM(n);
    end
        
    % Complete missing fields with defaults 
    % Don't report defaults yet, first calculate dependent fields according to new defaults...
    def = setdiff(fieldnames(DEF.(C.model)),fieldnames(C));
    C = completestruct(C,DEF.(C.model));
    
    % ... e.g. reference conditions according to IEC 61215 standard
    if C.isCPV, Gnoc = 900; VWnoc = 2; Tnoc = 20;   % CPV modules
    else, Gnoc = 800; VWnoc = 1; Tnoc = 20;         % non-concentrating
    end
    
    % ... or NOCT-based heat transfer coefficients
    if contains(C.model,{'Uval','Mattei'}) && isfield(C,'NOCT') && ismember('Uconst',def)   
        C.Uconst = C.absort*Gnoc/(C.NOCT - Tnoc);
        if C.Uwind > 0
            C.Uconst = C.Uconst - C.Uwind*VWnoc;
        end
        def = setdiff(def,{'Uconst'});
    else
        def = setdiff(def,{'NOCT'});
    end
    
    if ~strcmp(C.model,'SAPM'), def = setdiff(def,{'type'}); end
    if ~isempty(def)
        warning('getcelltempfcn:default',...
                'Using default cell-temperature parameters:\n%s\n',...
                dispnested(rmfield(C,setdiff(fieldnames(C),def))));
    end

    switch C.model
        case {'Uval'}
            if ~isfield(C,'eff'), C.eff = 0; end % NOCT
            if isnumeric(C.eff)
                ctf = @(g,t,v) celltempUvalues(g,t,v,C,C.eff);
            else
            % Create interpolant based on default grid. FUTURE: is this still necessary??
                ta = getSimOption('Tagrid');
                vw = getSimOption('vwgrid');
                ge = getSimOption('Gegrid');

                F = @(g,t,v) celltempUvalues(g,t,v,C,C.eff,TOL_TEMP/4);
                [~,~,F] = fitinterpolant(F,{ge,ta,vw},TOL_TEMP);
                ctf = @(g,t,v) F(g,t,v);
            end
            ctf0 = @(g,t,v) celltempUvalues(g,t,v,C,0);
        case 'SAPM'
            ctf = @(Ge,Ta,Vw) Ta + Ge.*exp(C.a_wind+C.b_wind*Vw) + C.delT*(Ge./1000);
            ctf0 = ctf;
        case 'Mattei'
            ctf = @(Ge,Ta,Vw) celltemp_mattei(Ge,Ta,Vw,C.Uconst,C.Uwind,C.absort,C.eta_ref,C.beta,C.gamma);
            ctf0 = @(Ge,Ta,Vw) celltemp_mattei(Ge,Ta,Vw,C.Uconst,C.Uwind,C.absort,0,0,0);
    end

    NOCT = ctf0(Gnoc,Tnoc,VWnoc);
    if isfield(C,'NOCT')
        if abs(NOCT - C.NOCT) > RELTOL/0.004
           warning('Nominal NOCT (%0.2f°C) does not quite fit cell-temp. model (%0.2f°C)',C.NOCT,NOCT);
        end
    else
        C.NOCT = NOCT;
    end
    
    C.info = celltempmodelstr(C);

    function S = checkmodelname(S,varargin)
        if isfield(S,'TempModel'), S.model = S.TempModel; end % compatibility with ODM
        if ~isfield(S,'model'), S.model = '(missing)'; end
        [S.model,found] = parselist(S.model,ALIAS,'cell-temperature model','',varargin{:});
        if ~found, S = rmfield(S,'model'); end
    end
    function S = parsemdl(S,model,src)
                
        % CHECKS{j,:} = {subset of KFIELDS}, @(x) field is valid, -quiet
        CHECKS = {{'model','type','info'},@(x) ischar(x) && ~isempty(x),1;
                  {'isCPV'},@(x) isscalar(x) && (isequal(x,0) || isequal(x,1)),0;
                  {'eff'}, @(x) (isa(x,'function_handle')) || ...
                    (isnumeric(x) && isscalar(x) && isreal(x) && (x >= 0) && (x <= 1)),0};
        CHECKS(end+1,:) = {setdiff(KFIELDS,cat(2,CHECKS{:,1})),...
                           @(x) isnumeric(x) && isscalar(x) && isreal(x),0};

        fld = [{'model';'NOCT'}; fieldnames(DEF.(model))];
        S.model = model;
        ignored = setdiff(fieldnames(S),fld);
        S = rmfield(S,ignored);
        ignored = intersect(ignored,KFIELDS);
        if ~isempty(ignored)
            warning('getcelltempfcn:irrelevant','Cell-temp. model = %s, ignoring %s from %s',...
                S.model,shortliststr(ignored,'irrelevant field','quotes',''''),src);
        end
        fld = fieldnames(S);
        invalid = false(numel(fld),1);
        for j = 1:size(CHECKS,1)
            idx = ismember(fld,CHECKS{j,1});
            if ~any(idx), continue; end
            bad = ~cellfun(@(f) CHECKS{j,2}(S.(f)),fld(idx));
            if CHECKS{j,3}
                idx(idx) = bad;
                S = rmfield(S,fld(idx));
            else
                invalid(idx) = bad;
            end
        end
        if any(invalid)
            warning('getcelltempfcn:invalid','Cell-temp. model = %s, ignoring %s from %s',...
                S.model,shortliststr(fld(invalid),'invalid field','quotes',''''),src);
        end
    end
end
function [ok,P] = testefficiency(eff)
    ok = isa(eff,'function_handle');
    if ~ok, P = struct(); return; end
    [g,tc] = meshgrid(10:400:1610,-20:20:60);
    eta = eff(g,tc);
    ok = all(eta(:) >= 0 & eta(:) <= 1);
    
    if nargout > 1
        P.eta_ref = eff(1000,25);
        b = [tc(:) - 25,log10(g(:)/1000)] \ (eta(:)./P.eta_ref - 1);
        [P.beta,P.gamma] = deal(-b(1),b(2));
    end
end
function s = celltempmodelstr(C)
% Take a CellTemp model structure with fields C.model {'Uval','SAPM'} and C.parameters (3-vector)
% and generate a user-feedback string describing the model. 

    switch C.model
        case {'Uval','NOCT'}
            s = sprintf('Faiman: Uc = %0.1f W/m²K, Uw = %0.1f W·s/m³K, a = %0.2f', ...
                C.Uconst,C.Uwind,C.absort);
        case 'SAPM'
            s = sprintf('SAPM: a = %0.1f ln(K·m²/W), b = %0.1f ln(K·m·s/W), dT = %0.1f °C', ...
                C.a_wind,C.b_wind,C.delT);
        case 'Mattei'
            s = sprintf(['Mattei: Uc = %0.1f W/m²K, Uw = %0.1f W·s/m³K, a = %0.2f, ', ...
                         'eta_ref = %0.1f%%, beta = %0.2f%%/K, gamma = %0.2f'],...
                C.Uconst,C.Uwind,C.absort,100*C.eta_ref,C.beta*100,C.gamma);
        otherwise
            keyboard();
    end

    if isfield(C,'NOCT'), s = sprintf('%s, NOCT = %0.1f °C',s,C.NOCT); end
    if isfield(C,'type'), s = sprintf('%s (%s)',C.type,s);
    elseif isfield(C,'isCPV') && C.isCPV, s = ['(CPV) ' s]; 
    end
end