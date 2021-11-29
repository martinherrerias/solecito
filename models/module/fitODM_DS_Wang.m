function P = fitODM_DS_Wang(Pmp,Vmp,Imp,Voc,Isc,Ki,Kv,Km)
% [Pmp,Vmp,Imp,Voc,Isc,c] = ONEDIODEMPP(PARAMS,TOL,MAXITER)  - Find the MPP of a one-diode model 
%   function defined by PARAMS, within precision TOL. Return also Voc, Isc, and c* as byproducts.
%
%  PARAMS: (array of) structure(s) of 5-8 ODM parameters (see ONEDIODEMODEL for details)
%  TOL: (incomplete) unified tolerance vector. See PARSETOLERANCE
%  MAXITER: optional limit of iterations for Newton-Raphson search.
%
%  [Pmp,Vmp,Imp,Voc,Isc,c] - size PARAMS arrays of MPP power, voltage, and current; open-circuit
%  voltage, short-circuit current, and curvature parameter c*, respectively.
%
%  (*) Parameter c = Vmp·(d²I/dV²)/(dI/dV) = -Vmp²/Imp·(d²I/dV²) @ MPP, can be used to get a 2nd 
%  order approximation of the IV curve in the vecinity of Vmp, namely:
%
%       I/Imp = a + b·exp(c·V/Vmp), with a = (1+1/c), and b = exp(-c)/c
%
%  c is useful to estimate the impact of small mismatch losses, proportional to (c+2)/c. See
%  L. L. Bucciarelli, “Power loss in photovoltaic arrays due to mismatch in cell characteristics” 
%  Solar Energy, vol. 23, no. 4, pp. 277–288, Jan. 1979, doi: 10.1016/0038-092X(79)90121-X.
%
% See also: ONEDIODEMODEL, TRANSLATEODM, PARSETOLERANCE


% Normalize irradiance by Reference, work only with non-zero values
    Ge = Ge./Mod.G_ref;
    if any(Ge < -eps | Ge > 2 | Tc < -50 | Tc > 100)
        warning('translateODM:GeTc','Working conditions outside expected boundaries');
    end
    Ge = max(Ge,0);
    
    k = 8.61733034e-5;         % Boltzmann's constant / unit charge (V/K)
    Mod.Tc_ref = Mod.Tc_ref+273.15;  % All temperatures to Kelvin
    Tc = Tc+273.15;
    Eg = Mod.Eg_ref*(1+Mod.muEg*(Tc - Mod.Tc_ref));
    n = Mod.nVth_ref/(k*Mod.Tc_ref*Mod.Ns);   % IdealityFactor·Njunctions

	if isfield(Mod,'muNdiode')
		n = n*(1+Mod.muNdiode*(Tc-Mod.Tc_ref));
	end
	[pars.nVth] = dealmat(n.*k.*Tc*Mod.Ns);
	 
    [pars.Iph] = dealmat(Mod.Iph_ref*Ge.*(1+Mod.muIsc*(Tc - Mod.Tc_ref)));
    [pars.Io] = dealmat(Mod.Io_ref*(Tc/Mod.Tc_ref).^3.*exp((Mod.Eg_ref./Mod.Tc_ref-Eg./Tc)./(n*k)));
    [pars.Rs] = deal(Mod.Rs);
	
	% Exponential Rsh correction
	if all(isfield(Mod,{'Rsh_0','Rsh_exp'}))
        [pars.Rsh] = dealmat(Mod.Rsh_ref+(Mod.Rsh_0-Mod.Rsh_ref)*exp(-Mod.Rsh_exp*Ge));
    else
        [pars.Rsh] = deal(Mod.Rsh_ref); %./Ge; 
	end
	
	% Pass parameters for recombination current, if available
	if all(isfield(Mod,{'di2mutau','Vbi'}))
		[pars.di2mutau] = deal(Mod.di2mutau);
		[pars.Vbi] = deal(Mod.Vbi);
	end
	
	% Pass reverse current characteristic factor, if available
	if isfield(Mod,'Brev')
		%[pars.Brev] = deal(Mod.Brev*Mod.Np/Mod.Ns^2);
        [pars.Brev] = deal(Mod.Brev/1000*Mod.Np/Mod.Ns^2); % mA/V² -> A/V²
	end


    % Get Voc,Isc with increased precission
    Voc = onediodemodel2(Params,0,tol);
    Isc = onediodemodel(Params,0,tol);

    % There is an actual maximum at Vmp < 0 for Voc < 0 and Isc < 0, but who cares
    if Isc <= toli(0) || Voc <= tolv(0)
        Vmp = 0; Imp = Isc; Pmp = 0; Voc = 0; Isc = 0; c = 0;
        return
    end
    
    % Debug
    % figure
    % plot(0:Voc/99:Voc,onediodemodel(Params,0:Voc/99:Voc));
    % hold on
    
    if ~isempty(Vmp_ref) && Vmp_ref > 0 && Vmp_ref < Voc
        Vmp = Vmp_ref;
    else
    % NOTE: 0.81·Voc is the median for most modules, but convergence is faster fromt the right
        Vmp = Voc*0.84; % Seed value 
    end
    
    [Imp,di,d2i] = onediodemodel(Params,Vmp,tol);
    for iter = 0:MaxIter 
        oldv = Vmp; %oldp = Pmp;
          
        if iter < 2
        % MPP for the approximation i/i0 = a - b·exp(cv/v_ref) [inspired by Bucciarelli, 1979]
        % NOTE: once LAMBERTWLOG is optimized/translated, it might be convenient to use this
        % approximation all the way. 
            c = Vmp*d2i/di;
            k = (1-d2i*Imp/di^2);
            Vmp = (lambertwlog(log(k)+c+1)-1)*Vmp/c;
        else
        % 2nd order Taylor expansion solution
            Vmp = Vmp-(Vmp*di+Imp)/(Vmp*d2i+2*di);
        end
        [Imp,di,d2i] = onediodemodel(Params,Vmp,tol);
        
        if abs(oldv-Vmp) <= tolv(Vmp) && abs((oldv-Vmp)*di) <= toli(Imp), break; end
    end
    if iter >= MaxIter
        % Vmp = NaN; Imp = NaN; Pmp = NaN;
        error('onediodeMPP:itermax','Maximum number of iterations reached');
    end
    Pmp = Vmp*Imp;
    c = Vmp*d2i/di;
end

function def = getdefaults(varargin)
% DEF = GETDEFAULTS(VARARGIN) - Parses NAME,VALUE list of options and returns a structure
% with fields {type, Eg_ref, muEg, nJunct, Vbi, Rsh_0_Rsh_ref, Rsh_exp, nDiode and di2mutau}.
% 
% Eg_ref = band_gap_ref (see https://duramat.github.io/pv-terms/)
% muEg = band_gap_temperature_coeff (id)
% nDiode = diode_factor (id)
%
% Default values are defined in ODMdefaults.mat
    %%
    [def,varargin] = getpairedoptions(varargin,{'type'},{[]});
   
    if isempty(def.type)
        warning('Using defaults for polycristalline silicon');
        def.type = 'pSi';
    end
    
    load('ODMdefaults.mat','materials');
    try def = materials(def.type,:);
    catch
        def = materials(contains(materials.tags,def.type),:);
        assert(size(def,1) == 1,'Unknown/ambiguous cell type');
    end
    type = def.Properties.RowNames{1};
    def = table2struct(def(1,2:end));
    
    def = getpairedoptions(varargin,def,'restchk');
    def.type = type;
end
