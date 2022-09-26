function Tc = celltempUvalues(Ge,Ta,vw,Mod,eff,Tol)
% Tc = CELLTEMPUVALUES(Ge,Ta,vw,Mod,eff,Tol)
%   Solves the Faiman Temperature Model: Tc = Ta + alpha·Ge·(1-eff)/(Uc + Uw·vw) for a Module with 
%   parameters [Uc,Uw,alpha] in Mod and efficiency eff(Ge,Tc), at conditions Ge,Ta,vw.
%
% Mod: 3-vector [Uc,Uw,alpha] or structure with fields 'Uwind', 'Uconst' and 'absort' 
%      Default values for Uwind and absort (if missing) are 0 W/(m²K·m/s) and 0.9.
% eff: function handle for efficiency: eff(Ge,Tc), or scalar, for eff(Ge,Tc) = e0
% Ge,Ta,vw: Effective irradiance [W/m²], ambient temperature [°C], and wind speed [m/s] arrays.
%           vw can be a scalar, that will be replicated for all elements of Ge, Ta
% NOTE: SI units are expected everywhere: W/m²,°C, m/s, W/m²K, W/(m²K·m/s). Use any other system at
% your own risk.
%
% References
%  D. Faiman, "Assessing the outdoor operating temperature of photovoltaic modules,"
%       Progress in Photovoltaics: Research and Applications, vol. 16, no. 4, pp. 307–315, 2008.
%
% See also: PVL_SAPMCELLTEMP, CELLTEMPUVALUES 

    narginchk(2,6);

    % errP/P = errT·(dP/dT)/P, for (dP/dT)/P = -0.004 1/K and errT = 100·tol -> |errP/P| = |0.4·tol|
    if nargin < 6 || isempty(Tol), Tol = getSimOption('RelTol')*100; end
    if nargin < 5 || isempty(eff), eff = 0.16; end
    if nargin < 4 || isempty(Mod), Mod = [29 0 0.9]; end
    if nargin < 3 || isempty(vw), vw = zeros(size(Ge)); end
    
    % If wind speed is a scalar, replicate for all Ee,Ta combinations
%     if isscalar(vw)&&numel(Ge)>1,vw = ones(size(Ge))*vw; end

    compatiblesize(Ge,Ta,vw);
%     if numel(Ge) ~= numel(Ta) || numel(Ge) ~= numel(vw)
%         error('celltempUvalues:arrsize','Ee, Ta, and/or vw have different sizes!');
%     end

    if isstruct(Mod)
        assert(isfield(Mod,'Uconst'),'celltempUvalues:missingUconst','Missing field ''Uconst'' in parameters');
        if ~isfield(Mod,'absort'), Mod.absort = 0.9; end
        if ~isfield(Mod,'Uwind'), Mod.Uwind = 0; end
        Mod = {Mod.Uconst,Mod.Uwind,Mod.absort};
    else
        assert(isvector(Mod) && numel(Mod) == 3, 'celltempUvalues:Mod','Expecting parameter structure or 3-vector');
        Mod = num2cell(Mod);
    end
    [uconst,uwind,alfa] = deal(Mod{:});

    if ~isa(eff,'function_handle')
    % No iterative solution required
       Tc = Ta + Ge.*(alfa-eff)./(uconst+uwind*vw);
    else
        % Get maximum temperature at eff = 0, i.e. cell is not reverse-biased
        Tmax = Ta + alfa*Ge./(uconst+uwind*vw); 

        % Solve by bisection method, note that Tmin = Ta @ eff = 1.
        %badeff = false;
        warning_resetter = naptime('bisection:closedintervals'); %#ok<NASGU>
        
        errf = @(t) Ta + Ge.*(alfa-effwrap(Ge,t))./(uconst+uwind*vw)-t;
        Tc = bisection(errf,Ta,Tmax,Tol);
    end
    
    function eta = effwrap(g,t)
        eta = eff(g,t);
        eta(g <= 0) = 0;
        assert(all(eta >= 0 & eta <= 1,'all'),'Bad efficiency function')
    end
end
