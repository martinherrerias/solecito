classdef BypassInterpolant
% Defines a Bypass Diodes's behavior from a set of ideal diode parameters. At construction it
% evaluates the model at different "ambient" (i.e. package < junction) temperatures and generates 
% an interpolation structure for fast IVcurve generation.
%
%	BypassInterpolant(pars,tc) contructs an object from diode parameters 'pars' at temperature
%							   values tc.
%	getIVpp(obj,tc)	returns the IV curve interpolated at the given temperature.
%
% See also: BYPASSDIODEFIT, BYPASSDIODE2, MODULEINTERPOLANT
	
properties (GetAccess = public, SetAccess = private)
    Pars	% Parameter structure {nDiode,Io,Rs,Imax,Vmax,Ir,...} see bypassdiode2 for details
    Tol     % tolerance used for construction [v,i,rel]
    Tc      % [1,N] vector of junction temperatures for reference curves [°C]
    Is      % [1,N] vector of temperature-adjusted leakage currents [A]
    Vmax    %
    Vmin    %
    Imax    %
    IVpp
end
methods
    % class constructor
    function obj = BypassInterpolant(S,Tc,Tol)
    % OBJ = BYPASSINTERPOLANT(PARS,TC,TOL) contructs a BYPASSINTERPOLANT object from 
    %   parameters PARS, evaluated at ambient (i.e. connection-box) temperatures TC.
    %
    % TC: vector of reference package temperatures [°C]
    % Tol: a 2-vector [TolV,TolI] for voltage and current tolerances, or a scalar to be
    %	multiplied by [Vmax,Imax] and used to generate PWL-approximations for the diode.
    
        if nargin < 2 || isempty(Tc), Tc = getSimOption('Tcgrid'); end
        if nargin < 3 || isempty(Tol) || any(Tol<=0), Tol = getSimOption('RelTol'); end
        
        if nargin < 1 || isempty(S), S = getSimOption('diode'); end
        
        if ischar(S)
        % Handle a mini-library of default diodes (for quick access)
            f = fileparts(mfilename('fullpath'));
            f = fullfile(f,[S '.diode']);
            if ~isempty(dir(f))
            % load existing file...
               D = load(f,'-mat'); D = D.obj;
            else
               D.Pars = 'foo'; 
            end
            S = bypassdiodefit(S,'forward');
            if isequal(S,D.Pars) && isequal(D.Tc,Tc)
            % ... if parameters match, use it
                obj = D;
            else
            % ... otherwise create interpolant, and save for next time
                obj = BypassInterpolant(S,Tc,Tol);
                save(f,'obj');
            end
            return;
        end
        
        % parse/complete parameters structure
        if any(~isfield(S,{'Io','nDiode','Rs','Tj0','XTI','Eg','Vmax','Imax','Ta0','U','Tjmax'}))
        
            S = bypassdiodefit(S);
        end
        assert(isfinite(S.Imax) && isfinite(S.Vmax),'Finite Vmax, Imax required');

        if isscalar(Tol), Tol = Tol*[S.Vfmax, S.Imax, 1]'; end
        Tol = parsetolerance(Tol);

        obj.Pars = S;
        obj.Tol = Tol;
        obj.Tc = Tc;
        
            % %DEBUG: Constant Tj approximations
            % nVth = @(tj) S.nDiode*8.61733034e-5*(tj+273.15);
            % Is = @(tj) bypassdiode(S,-Inf,tj);
            % muIs = @(tj) Is(tj)/(tj+273.15)*(S.Eg/nVth(tj)+1);
            % [tj,vv] = pwlapprox(@(tj) min(S.Vmax,-S.U./muIs(tj)),min(Tc),S.Tjmax,[0 Tol(1) Tol(3)]);
            % Tr = mdpwl(vv,tj);        
            % 
            % Tj = linspace(min(Tc),S.Tjmax,100);
            % [vv,tj] = meshgrid(linspace(-S.Vf,0,1000),Tj);
            % ii = bypassdiode(S,-vv,tj);
            % ta = tj +ii.*vv/S.U;
            % 
            % contourf(vv,ii,ta,Tc);
            % hold on;
            % plot(Tr.x,Is(Tr.y))
            % keyboard();

        Vmax = ones(size(Tc))*S.Vmax;
        obj.Is = -bypassdiode(S,-Vmax,Tc,Tol,'ambient');
        Vmin = bypassdiode2(S,S.Imax,Tc,Tol,'ambient');
        
        if any(isnan(obj.Is))           
            % Reverse voltage limited by thermal run-off (negative bias)
            for j = find(isnan(obj.Is))
                f = @(v) isnan(bypassdiode(S,v,Tc(j),0,'ambient'))-0.5;
                Vmax(j) = -bisection(f,0,-S.Vmax,0)-Tol(1);
                obj.Is(j) = -bypassdiode(S,-Vmax(j),Tc(j),0,'ambient');
            end 
        end
        obj.Vmax = Vmax;
        
        Imax = ones(size(obj.Is))*S.Imax;
        if any(isnan(Vmin))
        % Forward current limited by thermal run-off (forward bias)
            for j = find(isnan(Vmin))
                f = @(i) isnan(bypassdiode2(S,i,Tc(j),0,'ambient'))-0.5;
                Imax(j) = bisection(f,0,S.Imax,0)-Tol(2);
                Vmin(j) = bypassdiode2(S,Imax(j),Tc(j),0,'ambient');
            end
        end
        obj.Imax = Imax;
        obj.Vmin = Vmin;
        
        % remove restrictions, just to make sure there are no NaNs due to numerical precision
        S.Imax = Inf;
        S.Vmax = Inf;
        S.Tjmax = S.Tjmax + 1;

        % Get PWL approximations
        waitbarh = optwaitbar(0,'','Name','Creating/updating diode interpolant...'); 
        obj.IVpp = mdpwl.empty;
        N = numel(Tc);
        for k = N:-1:1            
            waitbarh.update(1-k/N,sprintf('Evaluating PWL approximation %d/%d',N-k,N));
            fun = @(v) bypassdiode(S,v,Tc(k),Tol,'ambient');
            [vv,ii] = pwlapprox(fun,-Vmax(k),Vmin(k),Tol);
            obj.IVpp(k) = mdpwl(-fliplr(vv),fliplr(ii),eps(0),[-Inf Vmax(k) -Inf Imax(k)]);
        end   
    end

    function IVpp = getIVpp(obj,t)
    % IVPP = OBJ.GETIVPP(T) - Return IV curve (MDPWL object) for junction temperature T (°C)

        X0 = obj.Pars.Vfmax;
        shift = @(p,dx) mdpwl(p.x+dx,p.y,Inf);
        
        W = interpmatrix(t(:),obj.Tc');
        assert(all(W(:) >= 0 & W(:) <= 1),'Extrapolation is not allowed')

        for j = numel(t):-1:1
            k = find(W(j,:) > 0);
            switch numel(k)
                case 1
                    IVpp(j) = copy(obj.IVpp(k));
                case 2
                    f = W(j,k(2));
                    pp = arrayfun(shift,obj.IVpp(k),[-X0 -X0]);
                    IVpp(j) = shift(interpolate(pp(1),pp(2),f),X0);
                    IVpp(j).clip([-Inf obj.Pars.Vmax -[1-f,f]*obj.Is(k)' [1-f,f]*obj.Imax(k)']);
                otherwise
                    error('You really should not be here, cry for help');
            end
        end
        IVpp = reshape(IVpp,size(t));
    end

    function newobj = scale(obj,Sv,Si)
    % Scale the IV curve bahavior described by obj by Sv in V, and by Si in I.
    % This can be used to represent arrays of obj elements (Sv,Si > 1 & integers),
    % sub-elements of obj (1/Sv,1/Si > 1 & integers), and/or reverse connections 
    % (Sv and/or Si == -1).

        newobj = obj;
        newobj.Ib = newobj.Ib*Si;
        newobj.Is = newobj.Is*Si;
        newobj.Vb = newobj.Vb*Sv;

        newobj.Pars.nDiode = obj.Pars.nDiode*Sv;
        newobj.Pars.Io = obj.Pars.Io*Si;
        newobj.Pars.Rs = obj.Pars.Rs*Sv/Si;
        newobj.Pars.Imax = obj.Pars.Imax*Si;
        newobj.Pars.Vmax = obj.Pars.Vmax*Sv;
        newobj.Pars.Ir = obj.Pars.Ir*Si;
        newobj.Pars.If = obj.Pars.If*Si;
        newobj.Pars.Vf = obj.Pars.Vf*Sv;
    end
    
    function plot(obj,tt)
    % OBJ.PLOT - generate forward-voltage vs current, and reverse leakage-current vs temperature
    %   plots for BYPASSINTERPOLANT object OBJ.
    
        ls = @(N,Xd,Xs) round((logspace(0,1,N)-1)/9*(Xs-Xd)+Xd);
        if nargin < 2, tt = ls(10,obj.Tc(end),obj.Tc(1)); end
                
        GUIfigure('Diode','Bypass Diode','2:0.9'); clf()
        set(gcf,'DefaultAxesColorOrder',parula(numel(tt)+1));
        subplot(1,2,1); hold on;
        title(obj.Pars.info);
        plot(obj.getIVpp(tt));       
        % plot(-obj.Vmin,obj.Imax,'r:');
        xlim([-round(max(obj.Vmin)*10+0.5)/10,0]);
        ylim([0 round(max(obj.Imax)+0.5)]);
        xlabel('(Forward) Voltage [V]');
        ylabel('(Forward) Current [A]');
        L = legend(num2str(tt(:),'%0.1f'));
        title(L,'T_{pckg} [°C]');
        set(L,'box','off');
        grid on
        
        subplot(1,2,2); hold on;
        plot(obj.getIVpp(tt));
        set(get(gca,'yaxis'),'exponent',-3);
        ylabel('(Reverse) Current [A]');
        
        % plot(obj.Vmax,-obj.Is*1000,'r:');
        xlim([0,round(max(obj.Vmax)/5+0.5)*5]);
        xlabel('(Reverse) Voltage [V]');
        grid on
    end
end

methods(Static)
    function obj = loadobj(obj)
    % obj = LOADOBJ(obj) - Generate STC parameters and interpolants on object load
    
        if ~isnestedfield(obj,'Is') || isempty(obj.Is)
        % Replace deprecated BYPASSINTERPOLANT objects
            fprintf('Updating BypassInterpolant object...\n')
            obj.Pars.Tjmax = 175;
            obj.Pars.RthJC = 4.0;
            obj.Pars.Tj = obj.Pars.Tj0;
            obj = BypassInterpolant(obj.Pars,[],obj.Tol);
            disp(obj);
        end
    end
end
end
