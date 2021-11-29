function pars = translateODM(ODM,Ge,Tc)
% PARS = TRANSLATEODM(ODM, Ge, Tc) - Calculate the adjusted One-Diode-Model parameters of model
%   ODM for effective irradiance(s) Ge and temperature(s) Tc'
%
% ODM - structure of parameters, typically returned by CHECKODM, with required numeric fields
%   {nVth_ref, Iph_ref, muIsc, Io_ref, Rsh_ref, Rs, Ns, Np, G_ref, Tc_ref, Eg_ref, muEg, muNdiode,
%   di2mutau,Vbi}, option keys {Rsh_model,reverse_model}, and, depending on these options, some of
%   the additional parameters {Rsh_0, Rsh_exp, Brev, Vbr, avalanche_a, avalanche_m}.
%
%   See CHECKODM and ONEDIODEMODEL for details.
%
% GE,TC - compatible size arrays of irradiance (units of ODM.G_ref) and cell temperature (°C)
%
% PARS - is a size(Ge) array of structures with fields {nVth,Iph,Io,Rsh,Rs} and,
%   depending on availability
%
% NOTE: Since translateODM is likely to be called over-and-over using the same ODM structure,  
%   this function leaves out most error checks! See CHECKODM.
%
% See also: ONEDIODEMODEL, ONEDIODEMODEL2, CHECKODM, MODULEINTERPOLANT

    narginchk(3,3);

    % Make sure inputs are of the same size
    [Ge,Tc] = compatiblesize(Ge,Tc);

    % Create output structure(s), set them to 0
    pars = repmat(struct('Iph',0,'Io',0,'Rs',0,'nVth',0,'Rsh',0),size(Ge));

    % Normalize irradiance by Reference, work only with non-zero values
    Ge = Ge./ODM.G_ref;
    if any(Ge < -eps | Ge > 2 | Tc < -50 | Tc > 100)
        warning('translateODM:GeTc','Working conditions outside expected boundaries');
    end
    Ge = max(Ge,0);
    
    k = 8.61733034e-5;  % Boltzmann's constant / unit charge (V/K)
    ODM.Tc_ref = ODM.Tc_ref+273.15;  % All temperatures to Kelvin
    Tc = Tc+273.15;
    Eg = ODM.Eg_ref*(1+ODM.muEg*(Tc - ODM.Tc_ref));
    n = ODM.nVth_ref/(k*ODM.Tc_ref*ODM.Ns);   % IdealityFactor·Njunctions

	if ODM.muNdiode ~= 0
		n = n*(1+ODM.muNdiode*(Tc-ODM.Tc_ref));
	end
	[pars.nVth] = dealmat(n.*k.*Tc*ODM.Ns);
	 
    [pars.Iph] = dealmat(ODM.Iph_ref*Ge.*(1+ODM.muIsc*(Tc - ODM.Tc_ref)));
    
    %                   I0=I0_ref.*((Tcell_K./Tref_K).^3).*exp((EgRef./(k.*Tref_K))-(E_g./(k.*Tcell_K)));
    [pars.Rs] = deal(ODM.Rs);
    
    [pars.di2mutau] = deal(ODM.di2mutau);
    [pars.Vbi] = deal(ODM.Vbi);
    
    switch lower(ODM.Io_model)
    case {'desoto','de soto','ideal'}
        [pars.Io] = dealmat(ODM.Io_ref*(Tc/ODM.Tc_ref).^3.*exp((ODM.Eg_ref./ODM.Tc_ref-Eg./Tc)./k));
    case {'pvsyst'}
        [pars.Io] = dealmat(ODM.Io_ref*(Tc/ODM.Tc_ref).^3.*exp((ODM.Eg_ref./ODM.Tc_ref-Eg./Tc)./(n*k)));
    otherwise
        error('Unrecognized reverse-saturation-current model');
    end
	
    switch ODM.Rsh_model
    case 'exponential'
        % PVsyst's Exponential Rsh correction
        [pars.Rsh] = dealmat(ODM.Rsh_base+(ODM.Rsh_0-ODM.Rsh_base)*exp(-ODM.Rsh_exp*Ge));
    case 'inverse'
        % De Soto, Klein & Beckman's inverse correction

        % NOTE: indiscriminate use of Rsh_ref/Ge leads to Rsh = Inf as Ge -> 0, which causes all
        %   sort of problems in reverse-bias models (as dI/dV @ Vj = 0 -> 0 ). Data from Ruschel
        %   et al. (2016), which supports inverse model over exponential correction, gets only 
        %   to 75 W/m² and registers values of Rsh/Rsh_ref of up to ~40. Setting a cap
        %   of 100 seems reasonable even for their most extreme case (Tandem), and will only
        %   affect curves with Ge < 10 W/m².

        [pars.Rsh] = dealmat(ODM.Rsh_ref.*min(100,1./Ge)); 
    end
    
    switch lower(ODM.reverse_model)
    case {'bishop','exponential'}
        [pars.avalanche_a] = deal(ODM.avalanche_a);
        [pars.avalanche_m] = deal(ODM.avalanche_m);
        [pars.Vbr] = deal(ODM.Vbr*ODM.Ns);
    case {'alonso','alonso-garcia','breakdown'}
    	if ~isfield(pars,'Vbi'), [pars.Vbi] = dealmat(ODM.Vbi); end
        [pars.Vbr] = deal(ODM.Vbr*ODM.Ns);
    case {'pvsyst','quadratic'}
        [pars.Brev] = deal(ODM.Brev*ODM.Np/ODM.Ns^2); % scale!
    case 'none'
    otherwise
        error('Unrecognized reverse-bias model');
    end
end

function varargout = dealmat(M)
% Deal the elements of M in varargout
    assert(nargout == numel(M), 'dealmat:narginNargoutMismatch','nargin-nargout mismatch');
    M = num2cell(M);
    varargout = M(:);
end
