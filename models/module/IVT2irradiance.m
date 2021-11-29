function Ge = IVT2irradiance(Mod,Ix,Vx,Tc,G_ref,Tol,MaxIter)
% Ge = IVT2irradiance(MOD,Ix,Vx,Tc,..) - Estimate effective incident irradiance on a PV element
%   modeled by One-Diode-Model MOD, given measured current(s) Ix, voltage(s) Vx, and cell-
%   temperature(s) Tc. That is, Find Ge such that:
%
%     Ix = ONEDIODEMODEL(P,Vx), for P = TRANSLATEODM(MOD,Ge,Tc)
%
% Ge = IVT2irradiance(MOD,Ix,Vx,fTc,G_ref) - Where fTc is a function handle, allows to consider the 
%   dependency of Tc on Ge. It iteratively solves for Ge, starting with seed G_ref (default Ix/Imp_ref)
%   and updating the guess for Tc = fTc(G) on every step. The choice of G_ref can be critical,
%   particularly for operational points close to Voc.
%
%   MOD - structure of extended one-diode-model parameters, as returned by CHECKODM.
%   Ix,Vx - equal sized arrays of current, voltage values.
%   Tc - size(Ix) array of cell-temperature values, or function handle @(G) Tc(G) that returns
%       such an array, given a guess for G.
%
% Ge = IVT2irradiance(MOD,Ix,Vx,Tc,G_ref,TOL,MAXITER) - Set tolerance/max. iteration options.
%
% See also: CHECKODM, ONEDIODEMODEL, TRANSLATEODM

    if nargin == 1 && isequal(Mod,'test'), IVT2irradiance_test(); return; end
    
    narginchk(4,7);
    
    G_MAX = 1600;
                
    if nargin < 7 || isempty(MaxIter), MaxIter = getSimOption('MaxIter'); end
    if nargin < 6 || isempty(Tol), Tol = getSimOption('RelTol'); end
    
    if nargin < 5 || isempty(G_ref), G_ref = Mod.G_ref*Ix/Mod.Imp_ref; end  % seed

    % Normalize irradiance by Reference, work only with non-zero values
    G_ref = min(max(G_ref,0),G_MAX)/Mod.G_ref;
    
    [Ix,Vx,G_ref] = compatiblesize(Ix,Vx,G_ref);
    
    if isa(Tc,'function_handle')
        fTc = Tc;
        Tc = fTc(G_ref*Mod.G_ref) + 273.15;
        last_Tc = Tc;
    else
        fTc = [];
        Tc = Tc+273.15;
    end
    compatiblesize(Vx,Tc);

    Vs = Vx + Mod.Rs*Ix;
    Vs = max(0,Vs); % return Isc for reverse-biased conditions

    % Recombination correction
    if all(isfield(Mod,{'di2mutau','Vbi'}))
        Fr = Mod.di2mutau./(Mod.Vbi - Vs);
    else
        Fr = 0;
    end

    k = 8.61733034e-5;  % Boltzmann's constant / unit charge (V/K)
    Mod.Tc_ref = Mod.Tc_ref+273.15;  % All temperatures to Kelvin
    
    for i = 1:MaxIter
    % Iterate to guess Tc, if required

        Eg = Mod.Eg_ref*(1+Mod.muEg*(Tc - Mod.Tc_ref));
        n = Mod.nVth_ref/(k*Mod.Tc_ref*Mod.Ns);   % IdealityFactor·Njunctions

        if isfield(Mod,'muNdiode')
            n = n*(1+Mod.muNdiode*(Tc-Mod.Tc_ref));
        end
        nVth = n.*k.*Tc*Mod.Ns;

        Iph_G_ref = Mod.Iph_ref.*(1+Mod.muIsc*(Tc - Mod.Tc_ref)); % Iph_G_ref
        Io = Mod.Io_ref*(Tc/Mod.Tc_ref).^3.*exp((Mod.Eg_ref./Mod.Tc_ref-Eg./Tc)./(n*k));

        if all(isfield(Mod,{'Rsh_0','Rsh_exp'}))
        % Exponential Rsh correction: Rsh = R_ref + (Rb - R_ref)e^(-bG)
        % Use 1st Taylor approx.: 1/Rsh ~ 1/Rsh·[1+b(1-R_ref/Rsh)(G-G_ref)] and iterate
        % Seems to converge quickly (2-3 iterations) for reasonable input values

            b = Mod.Rsh_exp;
            Rsh_ref = Mod.Rsh_ref;

            for j = 1:MaxIter
                Rsh = Mod.Rsh_ref+(Mod.Rsh_0-Rsh_ref)*exp(-b*G_ref);
                Rsh_slope = b*(1-Rsh_ref./Rsh); % Rsh·d(1/Rsh)/dG

                Ge = (Ix + Io.*(exp(Vs./nVth)-1) + Vs./Rsh.*(1-Rsh_slope.*G_ref))./...
                 (Iph_G_ref.*(1-Fr)-Rsh_slope.*Vs./Rsh);
             
                Ge = min(max(Ge,0),G_MAX/Mod.G_ref);

                if all(abs(Ge - G_ref),'all') < Tol, break; end
                G_ref = Ge;
            end
        else
        % Straightforward analytical solution
            Rsh = Mod.Rsh_ref;
            Ge = (Ix + Io.*(exp(Vs./nVth)-1) + Vs./Rsh)./(Iph_G_ref.*(1-Fr));
            Ge = min(max(Ge,0),G_MAX/Mod.G_ref);
        end
        
        if isempty(fTc), break;
        else
            Tc = fTc(Ge.*Mod.G_ref) + 273.15; % update cell-temp estimate
            if all(abs(Tc - last_Tc),'all') < Tol, break; end
            last_Tc = Tc;
            G_ref = Ge;
        end
    end

    Ge = Ge.*Mod.G_ref;
end

function IVT2irradiance_test()
    load(pickfile('*.ssp'),'-mat','celltemp','ModIVint');

    N = 1000;
    SIGMA.Tc = 3;
    SIGMA.V = ModIVint.Vmp_ref*0.02;
    SIGMA.I = ModIVint.Imp_ref*0.02;

    SIGMA.G = 100; % modeling uncertainty (seed value)

    G = rand(N,1)*1400;
    Ta = rand(N,1)*30;
    vw = rand(N,1)*5;

    Tc = celltemp(G,Ta,vw) + randn(N,1)*SIGMA.Tc;
    V = ModIVint.getVoc(G,Tc).*rand(N,1); % Let the system operate anywhere in [0,Voc]

    % I_pwl = arrayfun(@(g,t,v) ModIVint.getIVpp(g,t).val(v),G,Tc,V);
    I_ODM = arrayfun(@(p,v) onediodemodel(p,v),translateODM(ModIVint.source,G,Tc),V);

    V = V + randn(N,1)*SIGMA.V;
    I = I_ODM + randn(N,1)*SIGMA.I;

    % known Tc
    % Ge = IVT2irradiance(ModIVint.source,I,V,Tc);

    % Unknown Tc
    G_ref = min(max(0,G + randn(N,1)*SIGMA.G),1400);
    Ge = IVT2irradiance(ModIVint.source,I,V,@(g) celltemp(g,Ta,vw),G_ref);

    GUIfigure('IVT2Ge_test'); clf();
    scatter(G,Ge,1);
end

