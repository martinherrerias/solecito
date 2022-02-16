function ODM = checkODM(ODM)
% ODM = CHECKODM(ODM) - Evauate (and complete/adjust) a set of extended One-Diode-Model parameters 
% as described by (De Soto et al, 2006) and (Mermoud & Lejeune, 2010).
%
% PROVISIONAL: requires an almost complete set of parameters (See ODMSpec*.odm)
% FUTURE: merge an updated/improved version of COMPLETEODM. Include parameter-fitting & fine-
% tuning capabilities.
%
% ODM - structure with the following fields: 
%   [field - description (type/units) [default], [!] = required, no default available.]
%
%   BASIC DATA:
%     name - (string)
%     Nbpd - number of bypass-diodes [3]
%     width - (mm)
%     length - (mm)
%     material - {'c-Si';'a-Si';'CdTe';'CIGS';'CIS';'GaAs';'HIT';'GaInP/GaInAs/Ge'}
%                used for bandgap determination (Eg_ref,muEg), if missing. ['a-Si']
%     isCPV - Boolean [false]
%     area - Active PV area (m²) [!]
%     Ns - number of cells in series [!]
%     Np - number of cells in parallel [1]
%     nJunct - number of pn junctions [1]
%
%   STC SPECIFICATIONS (as calculated from model):
%     Voc_ref - Reference Open Circuit Voltage (V) [f(ODM)]
%     Isc_ref - Reference Short Circuit Current (A) [f(ODM)]
%     Imp_ref - Reference MPP current (A) [f(ODM)]
%     Vmp_ref - Reference MPP voltage (V) [f(ODM)]
%     Pmp_ref - Reference MPP power (W) [Imp_ref*Vmp_ref]
%     eta_ref - Reference MPP efficiency [Pmp_ref/(area·G_ref)]
% 
%     muPmp - MPP power dependence with temperature: d(Pmp/Pmp_ref)/dTc (1/K) [f(ODM)]
%     muVoc - Open Circuit Voltage dependence with temperature: d(Voc/Voc_ref)/dTc (1/K) [f(ODM)]
%
%   NOMINAL SPECIFICATIONS
%     nominal.x - for x in {Isc_ref, Imp_ref, Vmp_ref, eta_ref, Pmp_ref, Voc_ref, muPmp, muVoc} 
%     keeps original corresponding field (usually not exactly equal to its modelled counterpart). 
%
%     G_ref - Reference Irradiance (W) [1000]
%     Tc_ref - Reference Cell Temperature (°C) [25]
%
%   ONE-DIODE-MODEL BASE PARAMETERS
%     Rsh_ref - Reference Shunt Resistance (Ohm) [!]
%     Iph_ref - Reference photo-current (A) [!]
%     Io_ref - Reference saturation current (A) [!]
%     nVth_ref - Reference thermal voltage (V) [Ns·nDiode·nJunt·k·T/q]
%     nDiode - Diode quality factor [!]
%     Rs - Series resistnace (Ohm) [!]
% 
%   TEMPERATURE DERATE FACTORS
%     muIsc - Short Circuit Current dependence with temperature: d(Isc/Isc_ref)/dT (1/K) [!]
%     Eg_ref - Bandgap energy (at reference conditions) (eV) [f(material)]
%     muEg - Bandgap dependence with temperature: d(Eg/Eg_ref)/dTc (1/K) [f(material)]
%     muNdiode - nDiode temperature derate factor (1/K) [0]
%
%     Io_model - 'desoto' or 'pvsyst', the latter uses muEg = 0, and Eg_ref·nDiode for correction
%       of Io (see PVL_CALCPARAMS_PVSYST) The former allows any muEg, and uses Eg_ref directly
%       (see PVL_CALCPARAMS_DESOTO).
%
%   IRRADIANCE DEPENDENT PARAMETERS
%     Rsh_model - 'inverse' or 'exponential' forces the use of the inverse model of De Soto et al.
%       (1988) or the exponential correction of PVsyst. [default based on available fields]
%     Rsh_0 - Rsh(E=0), for exponential Rsh irradiance correction (Ohm) [4·Rsh_ref]
%     Rsh_exp - exponent for Rsh irradiance correction [5.5]
%
%   RECOMBINATION CURRENT 
%     Vbi - Built-in junction voltage, all cells (V) [0.9·Ns]
%     di2mutau - dimensionless recombination current term (di²/mu·tau) [0]
%
%   REVERSE-VOLTAGE CHARACTERISTIC
%     reverse_model - {'exponential','quadratic','breakdown','none'} forces the use of one of
%       the following reverse-bias models [default based on parameters available]
%
%     Bishop (1988) exponential model: Ia = Vj/Rsh·a·(1-Vj/Vbr)^(-m)  
%       avalanche_a - avalanche-fraction a [0.16]
%       avalanche_m - avalanche exponent m [2.27]
%       Vbr - breakdown voltage (per cell) [-18.5 V]
%
%     PVsyst quadratic approximation: Ia = Brev·Vj²
%       Brev - quadratic reverse-voltage coefficient (A/V²), single cell [1.2·Isc_ref]
%
%     Alonso-García & Ruiz (2006): i = i0(Vj)/(1-exp(Be·(1-sqrt((Vbi-Vbr)/(Vbi-V)))))
%         Vbr - breakdown voltage (per cell) [-18.5 V]
%         Vbi - Built-in junction voltage, all cells (V) [0.9·Ns]
%
%     No model
% 
% REFERENCES:
% De Soto, W., Klein, S.A., Beckman, W.A., 2006. Improvement and validation of a model for 
%   photovoltaic array performance. Solar Energy 80, 78–88.
% Mermoud, A., Lejeune, T., 2010. Performance assessment of a simulation model for PV modules 
%   of any available technology, in: 25th European PV Solar Energy Conference, Valencia.
% Sauer, K.J., Roessler, T., Hansen, C.W., 2015. Modeling the Irradiance and Temperature Dependence
%   of Photovoltaic Modules in PVsyst. IEEE Journal of Photovoltaics 5, 152–158.
% Merten, J., Asensi, J.M., Voz, C., Shah, A.V., Platz, R., Andreu, J., 1998. Improved 
%   equivalent circuit and analytical model for amorphous silicon solar cells and modules. 
%   Electron Devices, IEEE Transactions on 45, 423–429.
% Bishop, J.W., 1988. Computer simulation of the effects of electrical mismatches in 
%   photovoltaic cell interconnection circuits. Solar cells 25, 73–89.
% Alonso-García, M.C., Ruíz, J.M., 2006. Analysis and modelling the reverse characteristic of 
%   photovoltaic cells. Solar Energy Materials and Solar Cells 90, 1105–1120.
%
% See also: TRANSLATEODM, ONEDIODEMODEL, ONEDIODEMODEL2, ONEDIODEMPP

    % Required fields
    REQ = {'nDiode','Iph_ref','muIsc','Io_ref','Rsh_ref','Rs','Ns','area'};
    
    % Optional fields with unvariant defaults
    OPT = struct('Nbpd',3,'isCPV',false,'Tc_ref',25,'G_ref',1000,'material','c-Si','Np',1,...
        'nJunct',1,'muNdiode',0,'di2mutau',0,'Vbr',-18.5,'avalanche_a',0.16,'avalanche_m',2.27);
    
    % Defaults depend on other parameters, e.g. material, Ns, Np...
    DEP = {'nVth_ref','Eg_ref','muEg','Rsh_0','Rsh_exp','Rsh_base','Vbi','Brev'};

%     % Check Module for required fields
%     if isfield(ODM,'Isc0')
%         % try
%         % hmm...try updating from old naming convention(s) first
%             ODM = parsestruct(updateversion(ODM),REQ);
%         % end
%     end
    ODM = renameODMfields(ODM);
    parsestruct(ODM,REQ);

    % Pick model-variation defaults (partly) based on available parameters:
    original = fieldnames(ODM);
    fieldsarethere = @(varargin) isempty(setdiff(varargin(:),original));
    
    if ~fieldsarethere('Rsh_model')
        if fieldsarethere('Rsh_0','Rsh_exp'), OPT.Rsh_model = 'exponential';
        else, OPT.Rsh_model = 'inverse';
        end
    % else, use defaults, good luck
    end

    if fieldsarethere('Vbr','avalanche_a','avalanche_m'), OPT.reverse_model = 'exponential';
    elseif fieldsarethere('Brev'), OPT.reverse_model = 'quadratic';
    else
        if ~fieldsarethere('Vbi','Vbr')
            warning('checkODM:rev', 'No reverse-bias parameters, using default Vbi and/or Vbr')
        end
        OPT.reverse_model = 'breakdown';
    end

    OPT.Io_model = 'pvsyst'; % 'desoto'
    
    for f = fieldnames(ODM)'
       if isnumeric(ODM.(f{1})), ODM.(f{1}) = double(ODM.(f{1})); end  
    end
    
    % Fill in basic defaults and calculate dependent fields from others
    ODM = completestruct(ODM,OPT);
    ODM = getdependent(ODM);
    
    % Check constraints on ODM fields
    checkconstraints(ODM,REQ);
    
    % Re-parse model choices, and cleanup non-required parameters
    rmifdefault =  @(S,f) rmfield(S,setdiff(intersect(fieldnames(S),f),original));
    switch lower(ODM.Rsh_model)
    case {'pvsyst','exponential'}
    % PVsyst exponential correction model
        ODM.Rsh_model = 'exponential';
        ODM.Rsh_base = (ODM.Rsh_ref-ODM.Rsh_0*exp(-ODM.Rsh_exp))/(1-exp(-ODM.Rsh_exp));
        original{end+1} = 'Rsh_base'; % count implicit field as given
    case {'desoto','de soto','inverse'}
    % De Soto, Klein & Beckman: Rsh/Rsh_ref = G_ref/G
        ODM.Rsh_model = 'inverse';
        ODM = rmifdefault(ODM,{'Rsh_0','Rsh_exp','Rsh_base'});
    otherwise
        error('Unrecognized Shunt resistance model');
    end
    
    switch lower(ODM.reverse_model)
    case {'bishop','exponential'}
    % Bishop (1988) model
        ODM.reverse_model = 'exponential';
        ODM = rmifdefault(ODM,{'Brev'});
    case {'alonso','alonso-garcia','breakdown'}
    % Alonso-García & Ruiz (2006)
        ODM.reverse_model = 'breakdown';
        ODM = rmifdefault(ODM,{'avalanche_a','avalanche_m','Brev'});
    case {'pvsyst','quadratic'}
     % PVsyst
        ODM.reverse_model = 'quadratic';
        ODM = rmifdefault(ODM,{'avalanche_a','avalanche_m'});
    case 'none'
        ODM = rmifdefault(ODM,{'avalanche_a','avalanche_m','Brev','Vbr'});
    otherwise
        error('Unrecognized reverse-bias model');
    end

    new = rmfield(ODM,intersect(fieldnames(ODM),original));
    if ~isempty(fieldnames(new))
        struc2str = @(x) strjoin(strtrim(strsplit(evalc('disp(x)'),newline())),', ');
        warning('checkODM:assumption','ODM fields completed with defaults: %s',struc2str(new));
    end

    % Evaluate ODMtranslate over a grid of values, and check for negative parameters
    g = getSimOption('Gegrid'); g = linspace(g(1),g(end),8);
    t = getSimOption('Tcgrid'); t = linspace(t(1),t(end),8);
    [g,t] = meshgrid(g,t);
    P = translateODM(ODM,g,t);

    % Check constraints, this time on translated parameters
    checkconstraints(P,fieldnames(P));
    
    % For amorphous modules with recombination current term, check that Voc < Vbi
    if ODM.di2mutau > 0
        Voc = arrayfun(@(p) onediodemodel2(p,0),P);
        if ~all(Voc < ODM.Vbi)
            warning('checkODM:VbiVoc','Vbi < Voc for some operating conditions!')
        end
    end
   
    % Calculate STC parameters: 
    % {'Isc_ref','Imp_ref','Vmp_ref','eta_ref','Pmp_ref','Voc_ref','muPmp','muVoc'}
    try STC = getnominal(ODM);
    catch ERR
        error('ODM crashed at STC: \n\n%s',getReport(ERR));
    end
    nomfields = fieldnames(STC); 
    
    % Use calculated STC parameters, send any existing to substructure ODM.nominal
    for f = nomfields'
        if isfield(ODM,f{1}), ODM.nominal.(f{1}) = ODM.(f{1}); end
        ODM.(f{1}) = STC.(f{1});
    end
    
    % ... and check against nominal specs.
    T = [STC;cell2struct(repmat({NaN},numel(nomfields),1),nomfields)];
    T(2) = completestruct(ODM.nominal,T(2));
    T = struct2table(T,'RowNames',{'Model','Specs'});
    
    % Get error for each parameter...
    T{'Error (%)',:} = (1-T{'Model',:}./T{'Specs',:})*100;
    
    % Multiply all percentages x 100 (for display)
    f = {'muPmp','muVoc','eta_ref'};
    T{:,f} = 100*T{:,f};
    T{'Error (%)',f} = T{'Specs',f}-T{'Model',f};
    T.Properties.VariableNames(f) = cellfun(@(s) [s 'x100'],f,'unif',0); % 'X (%)' is not valid

    if any(abs(T{'Error (%)',:}) > 100*getSimOption('RelTol'))
        T{'Error (%)',:} = round(T{'Error (%)',:},4);

        warning('checkODM:nomchk',['ODM nominal specs are not fully consistent (|e|< %0.2f %%) '...
            'with model values. You might want to fine-tune your parameters to specs:\n\n%s'],...
            max(abs(T{'Error (%)',:})),evalc('disp(T)'));
    end
    
    function checkconstraints(S,requiredfields)
    % Check constraints on ODM fields
    
        TRANS = {'Iph','Io','Rs','nVth','Rsh'};  % returned by TRANSLATEODM
        STRINGS = {'name','material','TempModel','Rsh_model','reverse_model','Io_model'};
    
        % Define Groups that should be checked for constraints:
        allfields = [REQ,fieldnames(OPT)',DEP,TRANS,STRINGS];
        BOOLEAN = {'isCPV'};
        NUMERIC = setdiff(allfields,[BOOLEAN,STRINGS]);
        FINITE = setdiff(NUMERIC,{'Rsh_ref','Rsh'});
        NONPOSITIVE = {'Vbr'};
        NONNEGATIVE = setdiff(NUMERIC,{'muIsc','muNdiode','muEg','muPmp','muVoc'});    
        NONZEROS = {'Rsh_ref','area','nDiode','nJunct','Ns','nVth_ref','Eg_ref','Io','nVth','Rsh'};
        INTEGERS = {'Ns','Np','Nbpd','nJunct'};
    
        parseset = @(S,set,conditions) parsestruct(S,intersect(set,requiredfields),...
            'opt',setdiff(set,requiredfields),conditions{:});

        for fld = NONPOSITIVE
            if isfield(S,fld{1})
                for j = 1:numel(S), S(j).(fld{1}) = -S(j).(fld{1}); end 
            end
        end
                
        parseset(S,NUMERIC,{'numeric','real','scalar'});
        parseset(S,FINITE,{'numeric','finite'});
        parseset(S,NONNEGATIVE,{'numeric','nonnegative'});
        parseset(S,NONZEROS,{'numeric','nonzero'});
        parseset(S,INTEGERS,{'numeric','integer'});
        parseset(S,BOOLEAN,{'class',{'numeric','logical'},'binary'});
        parseset(S,STRINGS,{'class',{'char','string'}});
    end
end

function ODM = getdependent(ODM)
% Get dependent fields from others. Currently only {'nVth_ref','Eg_ref','muEg','Rsh_0','Vbi','Brev'}

    k = 8.61733034e-5;  % Boltzmann's constant (eV/K)
	
    % Fill in simple dependencies
    DEF = struct('Rsh_0',4*ODM.Rsh_ref,'Vbi',ODM.Ns*0.9,'Brev',1.2e-3*ODM.Iph_ref);
    ODM = completestruct(ODM,DEF);
    
    % Thermal voltage
    if ~isfield(ODM,'nVth_ref'), ODM.nVth_ref = ODM.nDiode*ODM.nJunct*k*(273.15+ODM.Tc_ref)*ODM.Ns; end
    
    % Bandgap lookup
    if ~isfield(ODM,'Eg_ref') || ~isfield(ODM,'muEg')
        
        % TODO: Unify material key names (later used for spectral correction EFFECTIVEIRRADIANCE)
        switch lower(ODM.material)
        case {'c-si','csi','monosi','simono','mono-c-si'}, key = 'c-Si';
        case {'p-si','psi','si-poly','multi-c-si'}, key = 'p-Si';
        case {'a-si','asi'}, key = 'a-Si';
        case {'cdte'}, key = 'CdTe';
        case {'cigs'}, key = 'CIGS';
        case {'cis'}, key = 'CIS';
        case {'gaas'}, key = 'GaAs';
        case {'hit','hit-si'}, key = 'HIT';
        case {'gainp','gainas','ge','gain'}, key = 'GaInP';
        otherwise
            key = 'c-Si';
            % ODM.material = [ODM.material ' (Eg = Si)'];
            warning('Unknown material: %s, using band-gap & spectral-response defaults for mono-Si',...
                ODM.material);
        end
        
        % Representative values, taken from PVLib v.1.1. PVL_CALCPARAMS_DESOTO:
        EGDB.types ={'c-Si';'p-Si';'a-Si';'CdTe';'CIGS';'CIS';'GaAs';'HIT';'GaInP'};
        EGDB.Eg_ref = [1.121;1.121;1.6;1.475;1.15;1.01;1.424;1.8;1.35];
        EGDB.muEg = [-0.0002677;-0.0002677;-0.0002677;-0.0003;0;-0.00011;-0.000433;-0.0002677;0];

        i = find(strcmp(EGDB.types,key));
        if ~isfield(ODM,'muEg'), ODM.muEg = EGDB.muEg(i); end    
        if ~isfield(ODM,'Eg_ref')
            ODM.Eg_ref = EGDB.Eg_ref(i); 
            ODM.Eg_ref = ODM.Eg_ref*(1 + ODM.muEg*(ODM.Tc_ref-25)); % bring 25°C value to Tc_ref
        end        
    end

    % FUTURE: fit missing parameters to nominal specs.
    %
    % if ~isfield(ODM,'Iph_ref'), ODM.Iph_ref = ODM.Isc_ref*(1+ODM.Rs/ODM.Rsh_ref); end    
    % if ~isfield(ODM,'nVth_ref')    
    %     if all(isfield(ODM,{'nDiode','Ns'}))
    %         ODM.nVth_ref = ODM.nDiode*ODM.nJunct*k*(273.15+ODM.Tc_ref)*ODM.Ns;
    %     else
    %         ODM.nVth_ref = (ODM.Voc_ref-ODM.Vmp_ref-ODM.Imp_ref*ODM.Rs)/...
    %             log((ODM.Isc_ref*(ODM.Rsh_ref+ODM.Rs)-ODM.Voc_ref)/...
    %                 ((ODM.Isc_ref-ODM.Imp_ref)*(ODM.Rs+ODM.Rsh_ref)-ODM.Vmp_ref));
    %         if ~isreal(ODM.nVth_ref)
    %             error('ODMchk:ComplexnVth',...
    %                   'Incompatible parameters, yield imaginary nVth');
    %         end
    %     end
    % end
    % 
    % if ~isfield(ODM,'Io_ref')
    %     ODM.Io_ref =-1/ODM.Rsh_ref*(ODM.Voc_ref-ODM.Isc_ref*(ODM.Rs+ODM.Rsh_ref))/...
    %         (exp(ODM.Voc_ref/ODM.nVth_ref)-1);
    % end
    % 
    % ODM.nDiode = ODM.nVth_ref/(k*(273.15+ODM.Tc_ref)*ODM.Ns*ODM.nJunct);
    % if ODM.nDiode < 1 || ODM.nDiode > 2 
    %     warning('ODMchk:DiodeNnotPhysical',...
    %             'Diode ideality factor (%0.2f) has no physical meaning',ODM.nDiode);
    % end
end

function M = getnominal(ODM)
    % Get MPP, Isc, Voc...
    odm_mpp = @(g,t) onediodeMPP(translateODM(ODM,g,t));
    [M.Pmp_ref,M.Vmp_ref,M.Imp_ref,M.Voc_ref,M.Isc_ref] = odm_mpp(ODM.G_ref,ODM.Tc_ref);
    M.eta_ref = M.Pmp_ref/(ODM.G_ref*ODM.area);
    
    % Temperature derate factors
    [M.muPmp,~,~,M.muVoc] = odm_mpp(ODM.G_ref,ODM.Tc_ref+1);
    M.muPmp = M.muPmp/M.Pmp_ref-1;
    M.muVoc = M.muVoc/M.Voc_ref-1;
end

 function M = renameODMfields(ODM)
% Rename fields from old file formats (backwards/CEC compatibility)

    KNOWN = {'name','Ns','Np','area','length','width','Nbpd','Io_model','nominal',...
    'Pmp_ref','Vmp_ref','Imp_ref','Isc_ref','Voc_ref','eta_ref','muPmp','muVoc',...
    'nDiode','Iph_ref','muIsc','Io_ref','Rs','nVth_ref','Vbi','di2mutau',...
    'Rsh_model','Rsh_ref','Rsh_0','Rsh_exp','Rsh_base',...
    'material','nJunct','isCPV','Eg_ref','muEg',...
    'reverse_model','Brev','Vbr','avalanche_a','avalanche_m',...        
    'NOCT','muNdiode','CEC_muIsc_adjust',...
    'TempModel','a_wind','b_wind','delT','Uconst','Uwind','absort',...
    'Tc_ref','G_ref','AM_ref','LibraryName','LibraryType'};

    REGEXP = {'([a-z])_?0$','$1_ref'; % X0 == X_0 == X_ref (*)
             '^([ivp])_mp','$1mp';  % I_mp == Imp, V_mp == Vmp, ..
             '^i_sc','Isc';         % I_sc == Isc
             '^v_oc','Voc';
             '^I0_','Io_';
             '^IL_','Iph_'};
     
    SIMPLE = {'rshbase','Rsh_0',1;
             'rshexp','Rsh_exp',1;
             'a_c','area',1;
             'gamma_r','muPmp',1./100;
             'gamma_ref','nDiode',1;
             'mugamma','muNdiode',1;
             'rs_ref','Rs',1;
             'ndiode_ref','nDiode',1;
             't_noct','NOCT',1;
             'source','material',1;
             'adjust','CEC_muIsc_adjust',1./100};
         
    warning('OFF','renamefields:overwrite');
    M = renamefields(ODM,REGEXP,'-regexprep','-ignorecase');
    M = renamefields(M,SIMPLE(:,1:2),'-ignorecase','scale',[SIMPLE{:,3}]);
    M = renamefields(M,[KNOWN;KNOWN]','-ignorecase');

    DEP = cell(0,3);
    if isfield(M,'Ns')
        DEP(end+1,:) = {'a_ref','nDiode',1./(8.61733034e-5 * 298.15 * double(M.Ns))};
    end
    if isfield(M,'Isc_ref')
        DEP(end+1,:) = {'alpha_sc','muIsc',1./M.Isc_ref};
    end
    if isfield(M,'Voc_ref')
        DEP(end+1,:) = {'beta_oc','muVoc',1./M.Voc_ref};
    end
    M = renamefields(M,DEP(:,1:2),'-ignorecase','scale',[DEP{:,3}]);
 end
