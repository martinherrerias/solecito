function [Tm,eta,pars] = celltemp_mattei(Ge,Ta,vw,varargin)
% [Tm,eta] = CELLTEMP_MATTEI(Ge,Ta,vw,[Uconst,Uwind,absort,eta_ref,beta,gamma])
% [Tm,eta] = CELLTEMP_MATTEI(Ge,Ta,vw,'Uconst',Uc,'Uwind',Uw,...)
% [Tm,eta,PARS] = CELLTEMP_MATTEI(Ge,Ta,vw,...)
%
%  Estimate module temperature and efficiency using the model of Mattei et al. (2006):
%
%     Upv = Uconst + Uwind·vw
%     eta = eta_ref·( 1 + beta·(Tm - 25°C) + gamma·log10(Ge/1000))
%     Tm = ( Upv · Ta + Ge·(absort - eta_ref·(1-beta·25°C) ) ) / ( Upv +  beta·eta_ref·Ge)
%
%  With:
%
%   Ge,Ta,vw: POA irradiance [W/m²], ambient temperature [°C], and wind speed [m/s] arrays
%
%   Uconst [ 24.1 W/m²K ] - Constant heat transfer coefficient 
%   Uwind [ 2.9 (W/m²K)/(m/s) ] - Wind dependent heat transfer coefficient 
%   absort [ 0.81 ]- Product of module absortivity and IR emmisivity (tau·alpha)
%   eta_ref [ 0.16 ] - Module efficiency at STC
%   beta [ 0.004 1/K] - (negative of) Pmp temperature derate factor (muPmp in ODM)
%   gamma [ 0 ] - Pmp irradiance log. derate factor (0-0.12 for silicon)
%
%  References:
%
% [1] Mattei, M., Notton, G., Cristofari, C., Muselli, M., Poggi, P., 2006. Calculation   
%   of the polycrystalline PV module temperature using a simple method of energy balance. 
%   Renewable Energy 31, 553–567. https://doi.org/10.1016/j.renene.2005.03.010
%
% See also: PVL_SAPMCELLTEMP, CELLTEMPUVALUES, GETCELLTEMPFCN

    narginchk(2,Inf);
    
    if nargin < 3 || isempty(vw)
    % complain, but do something, if there's no windspeed
        warning('Results will be terrible, you should get wind data!');
        vw = ones(size(Ge));
    end
    
    [Ge,Ta,vw] = compatiblesize(Ge,Ta,vw);

    % Define default options:
    pars.Uconst = 24.1;     % Intercept of heat transfer coefficient [W/m²K]
    pars.Uwind = 2.9;      	% Slope of heat transfer coefficient [(W/m²K)/(m/s)]
    pars.absort = 0.81;		% tau·alpha - Absortivity x IR Emissivity of module    
    pars.eta_ref = 0.16;    % efficiency at STC
    pars.beta = 0.004;      % Pmp temperature derate factor [1/K] (-muPmp in ODM)
    pars.gamma = 0.0;       % Pmp irradiance log. derate factor (typ. 0 - 0.12)

    pars = getpairedoptions(varargin,pars,'dealrest');

    % sanity checks
    parsestruct(pars,{'Uconst','Uwind','absort','eta_ref'},'-f','-r','-p','-s');
    validateattributes(pars.beta,'numeric',{'scalar','nonnegative','<',0.01});
    validateattributes(pars.gamma,'numeric',{'scalar','nonnegative','<',1});
    
    args = struct2cell(pars);
    [uconst,uwind,alpha,eta_ref,beta,gamma] = deal(args{:});
   
    T_STC = 25;
    
    Upv = uconst + uwind*vw;
    eta = eta_ref.*(1 + gamma*log10(Ge/1000) + beta*T_STC);
    
    Tm = (Upv.*Ta + Ge.*(alpha - eta))./(Upv - beta*eta_ref*Ge);
    if nargout > 1
        eta = eta - beta*eta_ref.*Tm;
    end
end
