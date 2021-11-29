function P = bypassdiodefit(S,varargin)
% S = BYPASSDIODEFIT(S) take (incomplete) specifications structure S and fit/complete all
%   parameters required for BYPASSDIODE and BYPASSDIODE2 to represent the model:
%
%       I = Is·(exp((V-Rs·I)/nVth) - 1),   for V >= -Vmax & I <= Imax
%
%     Where:
%         nVth = nDiode·k·T/q
%         Is = Io·(Tj/Tj0)^(XTI/nDiode)·exp(-Eg/(nDiode·k·Tj)·(1-Tj/Tj0))
%         Tj and Tj0 are the operational & nominal junction temperatures [K]
%         k/q is the Boltzmann constant / unit charge [8.61733034e-5 V/K]
%
%   REQUIRED fields:
%     If = N-Vector of forward current values [A] at Tj, Vf
%     Vf = N-Vector of forward current values [A] at Tj, If
%     Tj = N-Vector of junction temperatures [°C] for the points If,Vf
%     Imax = Max. forward current [A]
%     Vmax = Reverse breakdown voltage [V] (positive, i.e. -Vmin)
%     Tjmax = Max. junction temperature [°C]
%
%     At least one of Ta0, RthJC or U (see below)
%     At least two pairs Vf(Tj == Tj0), If(Tj == Tj0) are required to fit Rs, nDiode, otherwise
%       one (or both) must be provided.
%     At least two pairs Ir,Tj must be provided to fit XTI, Eg, otherwise... 
%
%   OPTIONAL fields:
%     Tj0 = [mode(Tj)] Nonminal junction temperature [°C]
%     Ir = reverse-bias (leakage) current [A] at -Vmax and Tj0 (scalar) or for each Tj (N-vector).
%     Ta0 = Package temperature [°C] @ Imax, Tjmax
%     U = Junction-to-case heat-transfer coefficient [W/K]
%     RthJC = Junction-to-case thermal resistance = 1/U
%
%   DEFAULT optional values (overwritten by fitted values, unless something goes wrong):
%     XTI = Saturation current temp. exponent (3.0 for pn-junction, 2.0 for Schottky-barrier)
%     Eg = Semiconductor energy bandgap [eV] (1.11 eV for Si)
%     Rs = Forward series resistance [Ohms] < Vf/If
%     nDiode = Diode ideality factor (theoretically 1-2)
%     Io = Saturation current at Tj0 [A]
%
%   OUTPUT structure S will have all of the fields listed above.
%
%   The algorithm for fitting parameters is the following:
%   
%       0. If Ir is provided, set Io = Ir(Tj0) (log-log interpolation)
%       1. Fit Rs, nDiode, and Io (if missing) from the subset of points If, Vf at Tj0
%       2. Fit temperature-derate parameters XTI, Eg based on all available If, Vf, Ir
%       3. Update Rs, nDiode to work with Eg, XTI, considering all points If, Vf, Ir
%       4. Go to 2. unless parameters have converged.
%
% S = BYPASSDIODEFIT(KEY,..) with KEY = 'pn', 'schottky' or 'smart' to use sample parameter sets.
%
% S = BYPASSDIODEFIT(..,'forward') - In practice, a diode has different saturation current behavior
%   at forward and reverse-bias conditions. A close fit on the reverse-bias region Is ~ f(Tj) will
%   likely result in poor temperature-derate performance on the forward-bias region. Using the
%   'forward' flag implies that steps 2 and 3 above ignore Ir constraints when fitting, ensuring a
%   better fit to the forward-bias region.
%
%   FUTURE: include separate Ior and nr parameters to model Is for V < 0
%   http://literature.cdn.keysight.com/litweb/pdf/ads15/ccnld/ccnld013.html
%
% S = BYPASSDIODEFIT(..,'plot') - Plot resulting Vf,If at junction temperatures Tj, and
%   reverse-leakage current vs temperature.
%
% See also: BYPASSDIODE2, CHECKODM, BYPASSINTERPOLANT

    if ischar(S)
        switch lower(S)
        case {'pn','silicon'}
        % Generic pn-junction diode
            P.info = 'Generic PN-junction diode';
            P.If = 10;
            P.Vf= 0.6;
            P.Ir = 2e-3;
            P.Tj = 125;
            P.Imax = 30;
            P.Vmax = 45;
            P.Tjmax = 175;
            P.Eg = 1.11;   % opt. parameters required for scalar If,Vf,Tj
            P.XTI = 3;     % ..
            P.Rs = 8e-3;   % ..
            P.RthJC = 4;
        case {'schottky','sbd'}
        % Microsemi Schottky-Diode SFDS1045Le3
            P.info = 'Schottky-Diode (Microsemi SFDS1045Le3)';
            P.Ir = [ 0.12 0.12  0.12 0.60   20  250  900]*1e-3;
            P.Tj = [   25   25    25   50  100  150  175];
            P.If = [    2    6    10   10   10   10   10];
            P.Vf = [0.343 0.39 0.425 0.40 0.35 0.31 0.29];
            P.Tj0 = 25;
            P.Imax = 10;
            P.Tjmax = 175;
            P.RthJC = 4.0; %[°C/W]
            P.Vmax = 45;
        case {'smart','fet'}
        % SM74611 Smart Bypass Diode (diode + FET)
            P.info = 'SM74611 Smart Bypass Diode';
            P.Ir = repmat([1e-4 4e-3 0.15 2],1,3)*1e-6;
            P.Tj = repmat([-40 25 85 125],1,3);
            P.If = [24 24 24 24 10 10 10 10 2 2 2 2];
            P.Vf = [56 75 105 133 22 30 47 67 7 10 20 47]*1e-3;
            P.Tj0 = 85;
            P.Imax = 15;
            P.Tjmax = 125;
            P.RthJC = 22.0; %[°C/W]
            P.Vmax = 28;
        otherwise
            error('Unknown model key');
        end
        S = P;
    end

    [~,TOL] = parsetolerance();
    MAX_ITER = 100;
    
    narginchk(1,3);
    opt = getflagoptions(varargin,{'plot','forward'},'restchk');
    P = S; % copy all specs, remember originals
    
    parsestruct(P,{'Imax','Vmax','Tjmax'},'-n','-r','-f','-s');
    parsestruct(P,{'If','Vf','Tj'},'-n','-r','-f','-e');
    parsestruct(P,{},'opt',{'Rs','nDiode','Eg','XTI','Tj0','U','RthJC','Ta0'},'-r','-f','-p','-s');

    pars = cellfun(@(x) P.(x),{'If','Vf','Tj','Imax','Tjmax'},'unif',0);
    [If,Vf,Tj,Imax,Tjmax] = deal(pars{:}); % Assign to individual variables
    If = If(:); Vf = Vf(:); Tj = Tj(:);
    
    Tj = Tj + 273.15;
    if ~isfield(P,'Tj0'), P.Tj0 = mode(P.Tj); end
    Tj0 = P.Tj0 + 273.15;
    Tjmax = Tjmax + 273.15;
    
    % Thermal Voltage = Temp [K] · Boltzmann constant k/q [V/K]
    Vth0 = Tj0*8.61733034e-5;
    Vth = Tj*8.61733034e-5;
    
    if isfield(P,'Ir')
        switch numel(P.Ir)
        case 1
           Ir = P.Ir;
           Io = Ir;
        case numel(Tj)
           Ir = P.Ir(:);
           [~,idx] = unique([Tj Ir],'rows');
           Io = exp(interp1(log(Tj(idx)/Tj0),log(Ir(idx)),0));
        otherwise
           error('Missing field Ir (scalar or %d-vector)',numel(Tj));
        end
    else
        idx = Tj == Tj0;
        assert(nnz(idx) > 2,'At least three points required with Tj == Tj0');
        b = [log(If(idx)).*Vth0,-Vth0*ones(nnz(idx),1),If(idx)]\Vf(idx);
        Io = exp(b(2)/b(1));
        Ir = [];
    end
    P.Io = Io;
    
    % Fit Rs, nDiode from a vector of points If,Vf @ Tj0
    idx = Tj == Tj0;
    A = [log(If(idx)./Io+1).*Vth0,If(idx)];
    nDiode = [];
    if rank(A) >= 2
        b = A\Vf(idx);
        if all(isreal(b) & b >= 0)      
           nDiode = b(1);
           Rs = b(2);
        end
    end
    if isempty(nDiode)
        if isfield(P,'nDiode')
            nDiode = P.nDiode;
            if isfield(P,'Rs')            
                Rs = P.Rs;
            else
                Rs = mean((Vf(idx) - nDiode*Vth.*log(1+If(idx)/Io))./If(idx));
            end
        else
            if isfield(P,'Rs')
                Rs = P.Rs; 
            else
                Rs = 0.5*min(Vf./If); 
            end
            nDiode = mean(((Vf(idx)-If(idx)*Rs)./Vth(idx))./log(1+If(idx)/Io));
        end
        P.nDiode = nDiode; P.Rs = Rs;
    end

    A = [log(Tj/Tj0),-(Tj0./Tj-1)/Vth0];
    if rank(A) < 2
        assert(all(isfield(P,{'XTI','Eg'})),...
            'Multiple junction-temperatures or explicit XTI, Eg required');
        Is = If./(exp((Vf-Rs*If)./(nDiode*Vth))-1);
    else
        prev = [nDiode Rs Inf Inf];
        for j = 1:MAX_ITER

            % Fit Eg/n and XTI/n from Is (both forward and reverse)
            Is = If./(exp((Vf-Rs*If)./(nDiode*Vth))-1);
            Is = max(Io,Is);
            if numel(Ir) > 1 && ~opt.forward
                b = [A;A]\log([Ir;Is]/Io);
            else 
                b = A\log(Is/Io);
            end    
            if all(isreal(b) & isfinite(b))      
               XTI_n = b(1);
                Eg_n = b(2);
            end

            % Correct Rs, nDiode using the new model for Is(Tj)
            Ism = Io.*(Tj/Tj0).^XTI_n.*exp(-Eg_n./Vth.*(1-Tj/Tj0));
            b = [log(If./Ism+1).*Vth,If]\Vf;
            if all(isreal(b) & b >= 0)      
                nDiode = b(1);
                Rs = b(2);
            end

            if all(abs([nDiode Rs Eg_n XTI_n] - prev) < TOL(prev)), break; end
            prev = [nDiode Rs Eg_n XTI_n];
            if j == MAX_ITER, error('Maximum number of iterations reached'); end
        end
        P.Rs = Rs; P.nDiode = nDiode; P.Eg = Eg_n*nDiode; P.XTI = XTI_n*nDiode;
    end

    P.Vfmax = bypassdiode2(P,P.Imax,-40);
    
    % Parse/Complete thermal-response parameters
    if isfield(P,'U')
        P.RthJC = 1/P.U;
        P.Ta0 = Tjmax -273.15 - (Imax*P.Vfmax)/P.U;
    elseif isfield(P,'RthJC')
        P.U = 1/P.RthJC;
        P.Ta0 = Tjmax -273.15 - (Imax*P.Vfmax).*P.RthJC;
    elseif isfield(P,'Ta0')
        assert(Tjmax - 273.15 > P.Ta0,'Max. junction temperature below ambient temperature');
        P.U = (Imax*P.Vfmax)/(Tjmax -273.15 - P.Ta0); 
        P.RthJC = 1/P.U;
    else
        error('U, RthJC or Ta0 are required to model thermal behavior of the diode')
    end
    
    % Issue warnings for parameters that have been overwritten
    overwritten = fieldnames(S);
    overwritten = overwritten(cellfun(@(f) ~isequal(P.(f),S.(f)),overwritten));
    if ~isempty(overwritten)
      warning('Replacing existing %s with fitted values',shortliststr(overwritten,'field')); 
    end

    if opt.plot
        [tj,~,ia] = unique(Tj);
        tj = tj - 273.15;
        ta = (25:25:Tjmax-273.15)';
        tj = [tj; ta(~ismember(round(ta/20),round(tj/20)))];

        GUIfigure('Diode','Diode','2:0.9'); clf();
        set(gcf,'DefaultAxesColorOrder',parula(numel(tj)+1));
        subplot(1,2,1); hold on;
        ii = logspace(-3,0,101)*Imax;
        arrayfun(@(t) plot(bypassdiode2(P,ii,t),ii),tj);
        scatter(Vf,If,20,ia);
        %xlim([-round(max(obj.Vmin)*10+0.5)/10,0]);
        %ylim([0 round(max(obj.Imax)+0.5)]);
        xlabel('(Forward) Voltage [V]');
        ylabel('(Forward) Current [A]');
        L = legend(num2str(tj(:),'%0.1f'),'Location','northwest');
        title(L,'T_j [°C]');
        set(L,'box','off');
        grid on

        subplot(1,2,2); hold on;
        tt = linspace(min(tj),max(tj),51);
        plot(tt,-bypassdiode(P,-Inf,tt)*1000,'r:');
        if numel(Ir) == numel(Tj)
            scatter(tj(ia),Ir*1000,20,ia,'Marker','x');
        end
        scatter(tj(ia),Is*1000,20,ia);
        set(gca,'yscale','log');
        xlabel('Junction Temperature [°C]');
        ylabel('Reverse Current [mA]');
        grid on    
    end
end

