function [Ix,Tj] = bypassdiode(S,Vx,varargin)
% Vx = BYPASSDIODE(PARS,Ix,TJ) - Return the forward voltage of a bypass diode defined by PARS
%   at forward voltage(s) Vx and junction temperature Tj. Solving the implicit equation: 
%
%     Ix = Is·(exp((Vx-Rs·Ix)/nVth) - 1),   for Vx >= -Vmax & i(Vx) <= Imax
%     Ix = NaN                              for Vx < -Vmax (diode breakdown)
%     Ix = Imax                             for i(Vx) > Imax (thermal failure)
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
%   Vx: array of voltage values [V]
%
% [Ix,Tj] = BYPASSDIODE(PARS,Vx,TA,TOL,'ambient') - Return the forward voltage of a bypass diode 
%   defined by PARS at forward voltage(s) Vx and AMBIENT temperature Ta. Instead of assuming a
%   known junction temperature, BYPASSDIODE(..,'-ambient') considers a constant junction-to-
%   environment heat transfer coefficient U and iteratively solves the system:
%
%     Tj = Ta[°C] + (Ix·Vx)/U + 273 K 
%     Ix = BYPASSDIODE(PARS,Vx,Tj)
%
%   PARS requires in this case additional parameters:
%
%     Tjmax = Max. junction temperature [°C]
%     Ta0 = Ambient temperature [°C] (for measurement of Imax, Tjmax)
%     U = Junction-to-ambient heat-transfer coefficient [W/K], U = (Imax·Vfmax)/(Tjmax - 25°C)
%     Vfmax = (forward) saturation voltage [V], Vfmax = bypassdiode2(PARS,Imax,Tjmax)
%
%   Vx: array of voltage values [V]
%   TA: ambient temperature [°C], if not specified, will be set to Ta0
%
%   As additional output you get TJ, the calculated junction temperature [°C] for each voltage.
%
%   NOTE: the system of equations above has two solutions for 0 > V > Vr, and no solutions for V
%   below Vr (diode burns due to thermal runaway). The runaway point occurs when the rate of
%   heat production excedes the rate of head dissipation, i.e.:
%
%       @ T = Tr, V = Vr:  d/dT[(T-Ta)U - IV] = 0  ->  U = V·dI(T,V)/dT
%
%   The runaway junction temperature Tr is thus an implicit funciton of voltage, and the search
%   space for the solution can be reduced to Ta < Tj < Tr(V).
%
% See also: BYPASSDIODE2, BYPASSDIODEFIT, BYPASSINTERPOLANT

    narginchk(2,5)

    if nargin > 2 && ischar(varargin{end})
        assert(isequal(lower(varargin{end}),'ambient'),'Uknown argument/flag');
        varargin = varargin(1:end-1);
        isambient = true;
    else
        isambient = false;
    end
    
    % parse common fields
    reqfields = {'Io','nDiode','Rs','Tj0','XTI','Eg'};
    parsestruct(S,reqfields(1:end-2),'-n','-r','-s','-p','-f');
    parsestruct(S,{'Vmax','Imax'},'-n','-r','-s','-p'); % let Vmax, Imax be Inf

    pars = cellfun(@(x) S.(x),reqfields,'unif',0);
    [Io,nDiode,Rs,Tj0,XTI,Eg] = deal(pars{:});
    Tj0 = Tj0 + 273.15;
    
    if ~isambient
    % Known junction temperature
        if  isempty(varargin) || isempty(varargin{1}), varargin{1} = Tj0; end
        
        Tj = varargin{1} + 273.15;
        Ix = diode(Tj,Vx);
        
        return;
    end
    
    % Unknown junction temperature: iterative solution:
    
    % parse addidional fields
    parsestruct(S,{'Ta0','U','Tjmax'},'-r','-s','-p','-f');    
    pars = cellfun(@(x) S.(x),{'Ta0','U','Tjmax'},'unif',0);
    [Ta0,U,Tjmax] = deal(pars{:});
    
    if isempty(varargin) || isempty(varargin{1}), varargin{1} = Ta0; end
    
    Ta = varargin{1} + 273.15;
    Tjmax = Tjmax + 273.15;

    % Initialize result arrays and unify sizes
    Ix = NaN(size(Ta+Vx));
    Tjmax = Tjmax.*ones(size(Ix));
    Vx = Vx.*ones(size(Ix));
    Ta = Ta.*ones(size(Ix));   
    
    warning_resetter = naptime({'bisection:closedintervals','bisection:notbounded'}); %#ok<NASGU>
    
    Tmin = Ta;
    Tmax = Tjmax +1; % give it some slack... (*)
    
%     % Reverse-bias: Find solution below thermal-runaway junction-temperature
%     f = Vx < (Tjmax-Ta)*U/diode(Tjmax,-Inf);
%     if any(f)
%         e = @(tj) Vx(f)-S.U./dIdT(tj,Vx(f));
%         Tmax(f) = bisection(e,Ta(f),Tjmax(f),0);
%     end
    
    % Find the thermal-runaway junction-temperature
    Q = @(t) (t-Ta)*U - diode(t,Vx).*Vx;
    dq = @(tj) Vx.*dIdT(tj,Vx)-S.U;
    [Tr,bound] = bisection(dq,Ta,Tmax,0);
    if any(bound)
        % Reverse-bias: Find solution below thermal-runaway junction-temperature
        f = bound & Vx < 0 & Q(Ta).*Q(Tr) <= 0;
        Tmax(f) = Tr(f);
%         f = ~f & Q(Tr).*Q(Tmax) <= 0;
%         Tmin(f) = Tr(f);
        
        % Forward-bias: Skip the trivial solution (Tj ~ Ta)
        f = bound & Vx > 0 & Q(Tr).*Q(Tmax) <= 0;
        Tmin(f) = Tr(f);
%         f = ~f & Q(Ta).*Q(Tr) <= 0;
%         Tmax(f) = Tr(f);
    end
    
    [Tj,bound] = bisection(Q,Tmin,Tmax,0);

%     % FUTURE: find threshold voltage using a fancier/faster method
%     if any(~bound)
%         for j = find(~bound)
%             f = @(t) diode(t,Vx(j)).*Vx(j)-(t-Ta(j))*U;
%             [Tmax(j),~,flag] = fminbnd(f,Ta(j),Tjmax(j));
%             if flag ~= 1
%                keyboard(); 
%             end
%         end
%         e = @(t) (t-Ta(~bound))*U - diode(t,Vx(~bound)).*Vx(~bound);
%     
%         [Tj(~bound),bound(~bound)] = bisection(e,Ta(~bound),Tmax(~bound),0);
%     end
    
    bound = bound & Tj <= Tjmax; % (*) ... and take it back

    Ix(bound) = diode(Tj(bound),Vx(bound));
    Tj = Tj - 273.15;

    function [Ix,dIdT] = diode(Tj,Vx)
        
        nVth = nDiode*8.61733034e-5*Tj; % n·Boltzmann constant k/q [V/K] · Temp [K]

        % Adjusted reverse saturation current
        %   Is = Io.*(Tj/Tj0).^(XTI/nDiode).*exp(-Eg./nVth.*(1-Tj/Tj0));
        logIs = log(Io)+(XTI/nDiode).*log(Tj/Tj0)-Eg./nVth.*(1-Tj/Tj0);
        Is = exp(logIs);
        
        if Rs == 0
        % Explicit solution
            Ix = Is.*(exp(Vx./nVth)-1);
        else
        % Solve using lambertw function:
            %   w = Is*Rs./nVth.*exp((Vx+Is*Rs)./nVth);
            %   Ix = nVth./Rs.*lambertw(w)-Is;
            logw = (Vx+Is*Rs)./nVth+logIs+log(Rs./nVth);
            Ix = nVth./Rs.*lambertwlog(logw)-Is;
        end
        
        % lambertwlog ensures precision within O(eps(W)) for W, but this doesn't imply sufficient
        % precission in Ix. Perform a few Newton-Raphson iterations directly over Ix, if necessary:
        % NOTE: use of exp(log(Io)+ ...) instead of Io·exp(...) helps reduce rounding error
            % g = @(Ix) exp(log(Is)+(Vx-Ix*Rs)/nVth)-Is - Ix;
            % df = @(Ix) -Rs/nVth*exp(log(Is)+(Vx-Ix*Rs)/nVth) - 1;
            % if any(abs(g(Ix)) > eps(Ix))
            %     Ix = quickfzero(g,df,Ix,@eps);
            % end
        
        % broke = Vx < -Vmax;  % Diode breakdown
        % Ix(broke) = NaN;   
        
        if nargout > 1
            dIdT = (Ix.*(Eg+nVth)-(Ix+Is).*(Vx-Rs*Ix))./(nVth+Rs*(Ix+Is))./Tj;
        end
    end

    function d = dIdT(Tj,Vx)
        [~,d] = diode(Tj,Vx);   
    end
end

