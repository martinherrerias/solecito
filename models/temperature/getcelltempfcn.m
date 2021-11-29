function [ctf,ctf0,C] = getcelltempfcn(S,eff)
% [ctf,ctf0,C] = GETCELLTEMPFCN(S,EFF) - Retrieve cell-temperature model parameters from structure
%   S. Check values, complete with defaults, and create cell-temperature-function/interpolants.
%
% [ctf,ctf0,C] = GETCELLTEMPFCN() - Interactive pick of SAPM temperature model, based on module
%   materials and mounting types.
%
% INPUT: structure S is expected to contain one of three sets of parameters:
%  
%     S.model = 'SAPM';  % Sandia Cell Temp. model [1]
%     S.a_wind = -3.47;  % Sandia Cell Temp. model coefficient 'a'
%     S.b_wind = -.0594; % Sandia Cell Temp. model coefficient 'b'
%     S.delT = 3.0;      % Sandia Cell Temp. model coefficient 'delta T' (K)
% 
%     S.model = 'Uval';  % Faiman quasi-physical model [2]
%     S.Uconst = 29.0;   % Heat transfer coefficient (wind independent) [W/m²K]
%     S.Uwind = 0.0;     % Heat transfer coefficient (wind dependent) [W·s/m³K]
%     S.absort = 0.9;    % optical absorptivity
% 
%     S.model = 'NOCT';  % Faiman model with Uconst adjusted to NOCT
%     S.NOCT = 44;       % Cell Temperature at NOC (°C)
%     S.Uwind = 0.0;   
%     S.absort = 0.9;
%
%   For the Faiman model (see below) a module efficiency function EFF(g,t) is required. If S is 
%   an ODM-parameter structure, however, the function can be generated internally as:
%
%       EFF = @(g,t) onediodeMPP(translateODM(ODM,g,t))./(ODM.area·g)
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
%
% See also: PVL_SAPMCELLTEMP, CELLTEMPUVALUES 

    RELTOL = getSimOption('RelTol');
    
    % For (dP/dT)/P ~ -0.004 1/K and dT = 100·tol -> |dP/P| = dT·(dP/dT)/P ~ |0.4·tol|
    TOL_TEMP = RELTOL*100;
    
    KFIELDS = {'model','Uconst','Uwind','absort','a_wind','b_wind','delT','isCPV','NOCT',...
              'TempModel','parameters'};
    
    % Copy only relevant fields
    if nargin < 1 || isempty(S), S = struct(); end
    C = rmfield(S,setdiff(fieldnames(S),KFIELDS));
    
    % Backwards compatibility: TempModel & parameters isntead of parameter names.
    if isfield(C,'TempModel'), C.model = C.TempModel; C = rmfield(C,'TempModel'); end 
    if isfield(C,'parameters')
       switch upper(C.model)
           case 'UVAL', f = {'Uconst','Uwind','absort'};
           case 'SAPM', f = {'a_wind','b_wind','delT'};
       end
       missing = ~cellfun(@(x) isfield(C,x),f);
       if any(missing)
           assert(numel(parameters) == 3,'');
           for j = find(missing)', C.(f{j}) = parameters(j); end
       end
       C = rmfield(C,'parameters');
    end
    
    % Load SAPM Cell-Temperature Model parameters, as default (just in case)
    [SAPM_PAR(1:6).model] = deal('SAPM');
    [SAPM_PAR(:).type] = deal('Glass/cell/glass - Open rack',...
                             'Glass/cell/glass - Close roof mount',...
                             'Glass/cell/polymer sheet - Open rack',...
                             'Glass/cell/polymer sheet - Insulated back',...
                             'Polymer/thin-film/steel - Open rack',...
                             '22X Linear Concentrator - Tracker');
    [SAPM_PAR(:).a_wind] = deal(-3.47,-2.98,-3.56,-2.81,-3.58,-3.23);
    [SAPM_PAR(:).b_wind] = deal(-.0594,-.0471,-.0750,-.0455,-.113,-.130);
    [SAPM_PAR(:).delT] = deal(3,1,3,0,3,13);
    
    % Load CellTemp structure from SimOptions, remember if it came from defaults
    [DEF_MDL, isdefault] = getSimOption('CellTemp');
    isdefault = all(cellfun(@(f) isdefault.(f),fieldnames(isdefault)));
    if nargin < 1 || isempty(C), C = DEF_MDL;
    else        
         % Resolve possible conflicts between CellTemp from S and from SimOptions
        fields = {'model','Uconst','Uwind','absort','a_wind','b_wind','delT'};
        fields = intersect(fields,fieldnames(DEF_MDL));
        if ~isdefault && ~all(cellfun(@(f) isequal(DEF_MDL.(f),C.(f)),fields))
            switch optquestdlg({['The Cell-Temperature model provided in AddSimOptions is '...
                    'inconsistent with that found in the Module-Structure:'];'';...
                    [' SimOptions: ' celltempmodelstr(DEF_MDL)];...
                    ['     Module: ' celltempmodelstr(C)];'';...
                    'Please choose one of the two.'},'Cell-Temperature Model',...
                    'SimOptions','Module','Module')
                case 'SimOptions', for f = fields(:)', C.(f{1}) = DEF_MDL.(f{1}); end
                case 'Module' % model = model;
                otherwise, error('Cancelled by user');
            end
        end
        isdefault = false;
    end
    assert(isfield(C,'model'),'Missing field ''model'' in (default?) structure.');

    % Get/check efficiency function (required for Faiman model)
    if ~strcmpi(C.model,'SAPM') 
        try 
            if nargin < 2 || isempty(eff)
                eff = @(g,t) onediodeMPP(translateODM(S,g,t),TOL_TEMP/4)./(S.area*g); 
            end
            [ge,tc] = meshgrid(1:800:1601,-20:40:60);
            eta = eff(ge,tc);
        catch, eta = -1; % crash below
        end
        assert(all(eta(:) >= 0 & eta(:) <= 1),'Missing or non-physical efficiency function');
        clear ge tc eta
    end

    % Reference conditions according to standards
    if isfield(C,'isCPV') && C.isCPV, Gnoc = 900; VWnoc = 2; Tnoc = 20;   % CPV modules
    else, C.isCPV = false; Gnoc = 800; VWnoc = 1; Tnoc = 20;              % non-concentrating
    end
    
    % interactive SAPM model pick
    if isdefault || (strcmpi(C.model,'SAPM') && ~isfield(C,'a_wind'))
        if runningfromUI()
            n = listdlg('PromptString','Choose module type and mounting','ListSize',[300,100],...
                    'ListString',{SAPM_PAR(:).type},...
                    'SelectionMode','single','Name','SAPM Cell-Temp. Model');
            assert(~isempty(n),'Invalid Module/mointing type!');
        else, n = 1;
        end
        for f = fieldnames(SAPM_PAR)', C.(f{1}) = SAPM_PAR(n).(f{1}); end
    end
    
    warnassump = @(msg,varargin) warning('ODMchk:assumption',['Assumption is made that ',msg],varargin{:});

    % Complete missing fields with defaults
    switch upper(C.model)
        case 'NOCT'
            try
                if ~isfield(C,'absort'), C.absort = 0.9; warnassump('absort = 0.9'); end
                if ~isfield(C,'Uwind'),C.Uwind = 0; warnassump('Uwind = 0'); end

                % eta_noc = 0! no electrical load is a condition of NOCT - IEC 61215
                C.Uconst = C.absort*Gnoc/(C.NOCT - Tnoc);
                if C.Uwind>0, C.Uconst = C.Uconst - C.Uwind*VWnoc; end
                C.model = 'Uval';
                warning('ODMchk:NOCT','Temp. model Uconst (%0.2f W/m²K) from NOCT',C.Uconst);
            catch ERR
                if strcmpi(DEF_MDL.model,'NOCT'), rethrow(ERR); end
                warning('ODMchk:NOCTfail','Failed to get Temp. model from NOCT, using defaults');
                [ctf,ctf0,DEF_MDL] = getcelltempfcn([],eff);
                for f = fieldnames(DEF_MDL)', C.(f{1}) = DEF_MDL.(f{1}); end
                return;
           end 
        case 'UVAL'
            C.model = 'Uval';
            if ~isfield(C,'absort'), C.absort = 0.9; warnassump('absort = 0.9'); end
            if ~isfield(C,'Uwind'),C.Uwind = 0; warnassump('Uwind = 0 W·s/m³K'); end
            if ~isfield(C,'Uconst'),C.Uconst = 29; warnassump('Uconst = 29 W/m²K'); end
        case 'SAPM'
            C.model = 'SAPM';
            if ~isfield(C,'b_wind')
                C.b_wind = SAPM_PAR(1).b_wind; 
                warnassump('b_wind = %0.2f ln(K·m·s/W)',C.b_wind); 
            end
            if ~isfield(C,'delT')
                C.delT = SAPM_PAR(1).delT; 
                warnassump('delT = %0.1f °C',C.delT); 
            end
        otherwise
            error('No valid Cell-Temp. model could be retrieved, either from SimOptions or ODM')
    end

    switch C.model
        case 'Uval'
            % Create interpolant based on default grid. FUTURE: is this still necessary??
            ta = getSimOption('Tagrid');
            vw = getSimOption('vwgrid');
            ge = getSimOption('Gegrid');

            F = @(g,t,v) celltempUvalues(g,t,v,C,eff,TOL_TEMP/4);
            [~,~,F] = fitinterpolant(F,{ge,ta,vw},TOL_TEMP);
            ctf = @(g,t,v) F(double(g),double(t),double(v));
            ctf0 = @(g,t,v) celltempUvalues(g,t,v,C,0);
        case 'SAPM'
            ctf = @(Ge,Ta,Vw) Ta + Ge.*exp(C.a_wind+C.b_wind*Vw) + C.delT*(Ge./1000);
            ctf0 = ctf;
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
end
    
function s = celltempmodelstr(C)
% Take a CellTemp model structure with fields C.model {'Uval','SAPM'} and C.parameters (3-vector)
% and generate a user-feedback string describing the model. 

    switch C.model
        case 'Uval'
            s = sprintf('Faiman: Uc = %0.1f W/m²K, Uw = %0.1f W·s/m³K, a = %0.2f', ...
                C.Uconst,C.Uwind,C.absort);
        case 'SAPM'
            s = sprintf('SAPM: a = %0.1f ln(K·m²/W), b = %0.1f ln(K·m·s/W), dT = %0.1f °C', ...
                C.a_wind,C.b_wind,C.delT);
        otherwise
            keyboard();
    end

    if isfield(C,'NOCT'), s = sprintf('%s, NOCT = %0.1f °C',s,C.NOCT); end
    if isfield(C,'type'), s = sprintf('%s (%s)',C.type,s);
    elseif isfield(C,'isCPV') && C.isCPV, s = ['(CPV) ' s]; 
    end
end