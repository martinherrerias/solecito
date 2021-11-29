function [Vx,Tj] = bypassdiode2(S,Ix,varargin)
% Vx = BYPASSDIODE2(PARS,Ix,TJ) - Return the forward voltage of a bypass diode defined by PARS
%   at forward current(s) Ix and junction temperature TJ. Using the equations: 
%
%     Vx = v, such that Ix = Is·(exp((v-Rs·Ix)/nVth) - 1),   for v(Ix) >= -Vmax & Ix <= Imax
%     Vx = -Vmax                                             for v(Ix) < -Vmax (diode breakdown)
%     Vx = NaN                                               for Ix > Imax (thermal failure)
%     
%   where:
%     nVth = nDiode·k·T/q, where T = Tj[°C] + 273 K 
%     Is = Io·(Tj/Tj0)^(XTI/nDiode)·exp(-Eg/(nDiode·k·Tj)·(1-Tj/Tj0))
%
%   PARS: parameter structure (see BYPASSDIODEFIT) with fields:
%
%     Io = Saturation current at Tj0 [A]
%     nDiode = Diode ideality factor (or emission coefficient)
%     Rs = Forward series resistance [Ohms]
%     Vmax = Reverse breakdown voltage [V]
%     Imax = Max forward current [A]
%     Tj0 = Nonminal junction temperature [°C]  
%     XTI = Saturation current temp. exponent (3.0 for pn-junction, 2.0 for Schottky-barrier)
%     Eg = Semiconductor bandgap [eV] (1.11 eV for Si)
%
%   Ix: array of current values [A]
%   Tj: junction temperature [°C], if not specified, will be assumed at Tj0
%
% [Vx,Tj] = BYPASSDIODE2(PARS,Ix,TA,TOL,'ambient') - Return the forward voltage of a bypass diode 
%   defined by PARS at forward current(s) Ix and AMBIENT temperature Ta. Instead of assuming a
%   known junction temperature, BYPASSDIODE2(..,'-ambient') considers a constant junction-to-
%   environment heat transfer coefficient U and iteratively approximates junction temperature as:
%
%     Tj = Ta[°C] + (Ix·Vx)/U + 273 K 
%     Vx = BYPASSDIODE2(PARS,Ix,Tj)
%
%   PARS requires in this case additional parameters:
%
%     Tjmax = Max. junction temperature [°C]
%     Ta0 = Ambient temperature [°C] (for measurement of Imax, Tjmax)
%     U = Junction-to-ambient heat-transfer coefficient [W/K], U = (Imax·Vfmax)/(Tjmax - 25°C)
%
%   Ix: array of current values [A]
%   TA: ambient temperature [°C], if not specified, will be set to Ta0
%
%   As additional output you get TJ, the calculated junction temperature [°C] for each voltage.
%
%   NOTE: The solution of the system might not be unique, at least for currents below Is(Ta),
%   the (negative) saturation current at junction temperature Ta. 
%   PROVISIONAL: there is no guarantee that the function will stick to the same solution branch, 
%   or that it will be consistent.
%
% See also: BYPASSDIODE, BYPASSDIODEFIT, BYPASSINTERPOLANT

    narginchk(2,5)

    if nargin > 2 && ischar(varargin{end})
        assert(isequal(lower(varargin{end}),'ambient'),'Uknown argument/flag');
        varargin = varargin(1:end-1);
        isambient = true;
    else
        isambient = false;
    end
    
    % parse common fields
    reqfields = {'Io','nDiode','Rs','Tj0','XTI','Eg','Vmax','Imax'};
    parsestruct(S,reqfields(1:end-2),'-n','-r','-s','-f');
    parsestruct(S,{'Vmax','Imax'},'-n','-r','-s','-p'); % let Vmax, Imax be Inf

    pars = cellfun(@(x) S.(x),reqfields,'unif',0);
    [Io,nDiode,Rs,Tj0,XTI,Eg,Vmax,Imax] = deal(pars{:});
    Tj0 = Tj0 + 273.15;
    
    if ~isambient
    % Known junction temperature:
        if  isempty(varargin) || isempty(varargin{1}), varargin{1} = Tj0; end
                
        Tj = varargin{1} + 273.15;
        Vx = diode(Tj,Ix);
        Vx(Ix > Imax) = NaN; % Diode (probably) burns

        return;
    end
    
    % Unknown junction temperature: iterative solution:
    
    % parse addidional fields
    parsestruct(S,{'Ta0','U','Tjmax'},'-r','-s','-p','-f');    
    pars = cellfun(@(x) S.(x),{'Ta0','U','Tjmax'},'unif',0);
    [Ta0,U,Tjmax] = deal(pars{:});

    if numel(varargin) < 2, tol = []; else, tol = varargin{2}; end
    tol = parsetolerance(tol);                          % [tv,ti,rel]
    % Tolerance vector adjusted to search function (*): 
    tol = parsetolerance([0,tol(1)*tol(2),tol(3)]);     % [0,tp,rel]

    if isempty(varargin) || isempty(varargin{1}), varargin{1} = Ta0; end
    
    Ta = varargin{1} + 273.15;
    Tjmax = Tjmax + 273.15;
    
    Vabsmax = max(0,Tjmax - Ta)*U./abs(Ix);
    burnt = abs(diode(Tjmax,Ix)) > Vabsmax;
    Vx = zeros(size(burnt));
    Tj = zeros(size(burnt));
    
    Tjmax = Tjmax.*ones(size(burnt)); Tjmax = Tjmax(~burnt);
    Ix = Ix.*ones(size(burnt)); Ix = Ix(~burnt);
    Ta = Ta.*ones(size(burnt)); Ta = Ta(~burnt);
    
    % (*): Search function f(T) = P_dissipated - P_produced = 0 (steady-state)
    f = @(t) (t-Ta)*U - Ix.*diode(t,Ix);
    Tj(~burnt) = bisection(f,Ta,Tjmax,tol);

    Vx(~burnt) = diode(Tj(~burnt),Ix);
    Vx(burnt) = NaN;     % Diode burns
    
    Tj(burnt) = NaN;
    Tj = Tj - 273.15;
    
    function Vx = diode(Tj,Ix)
        
        nVth = nDiode*8.61733034e-5*Tj; % n·Boltzmann constant k/q [V/K] · Temp [K]

        % Adjusted reverse saturation current
        % Is = Io.*(Tj/Tj0).^(XTI/nDiode).*exp(-Eg./nVth.*(1-Tj/Tj0));
        logIs = log(Io)+(XTI/nDiode).*log(Tj/Tj0)-Eg./nVth.*(1-Tj/Tj0);
        Is = exp(logIs);

        Vx = nVth.*log(Ix./Is+1)+Rs*Ix;

        broke = Ix < -Is | Vx < -Vmax;        % Set V = -Vmax @ Diode breakdown
        Vx(broke) = -Vmax;
    end
end
