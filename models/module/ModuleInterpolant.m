classdef ModuleInterpolant < matlab.mixin.Copyable
% Defines a PV-module's behavior from either a One-Diode-Model (ODM) parameter set, or directly 
% from a set of reference curves at diferent effective irradiance and cell temperature conditions.
	
properties
    name
    source  % model parameter structure (ODM / SAPM); empty for measured IVpps.
    info    % nominal parameter structure, free field names
    G0      % Reference Effective Irradiance (W/m²)
    Tc0     % Reference Cell Temperature (°C)
    area    % m²
    geom    % pvArea object, with dimensions and cell-connectivity
end
properties (GetAccess = public, SetAccess = private, Transient = true)
    Isc0    % Short circuit current [A], reference conditions
    Imp0    % Max. Power Point current [A], reference conditions
    Vmp0    % MPP voltage [V], reference conditions
    Voc0    % Open Circuit Voltage [A], reference conditions
end
properties (GetAccess = public, SetAccess = private)
    Gref	% vector of reference effective-irradiance values (W)
    Tref	% vector of reference cell temperatures (°C)
    Tol     % tolerance with which the object was generated (ODM) or precision of IVpps
    IVpps   % array of PWL curves (MDPWL objects) for ndgrid(Gref,Tref)
end
properties (Hidden = true, Access = private) % Transient = true
% Gridded interpolants f(g,t) accessible by get methods
    PmpInterp
    ImpInterp
    VocInterp
    IscInterp
    cInterp
end
methods(Access = protected)
    function cp = copyElement(obj)
    % Called by copy(obj), deep-copy of obj, including PVAREA and MDPWL objects

        cp = copyElement@matlab.mixin.Copyable(obj); % Shallow copy
        cp.geom = copy(obj.geom);
        cp.IVpps = arrayfun(@copy,obj.IVpps);
    end 
end
methods
    function obj = ModuleInterpolant(source,ge,tc,varargin)
    % OBJ = MODULEINTERPOLANT(IVPPS,GE,TC,PARAMS) contructs a MODULEINTERPOLANT object from 
    %   a (gridded) set of IV-curves IVPPS, evaluated at effective-irradiances GE and cell-
    %   temperatures TC. Optionally append nominal parameters in PARAMS.
    %
    % OBJ = MODULEINTERPOLANT(ODM,GE,TC,LIMS,TOL) first generates a set of IVPPS from the
    %   One-Diode-Model parameter structure ODM (see TRANSLATEODM and ONEDIODEMODEL), starting
    %   with (but not limited to) the set of conditions NDGRID(GE,TC). PWL approximations are 
    %   clipped to limits LIMS and calculated with tolerance TOL.
    %
    % Input:
    %   IVPPS - an Ng·Nt array of MDPWL curves for conditions ndgrid(GE,TC)
    %   GE - [Seed] vector of reference irradiances (Wh/m²)
    %   TC - [Seed] vector of reference cell-temperatures (°C)
    %   PARAMS - nominal parameter structure with free-field names. Specially handled 
    %       will be {'name','area','G0','Tc0','geom'}, all others will be moved to OBJ.info
    %   LIMS - 4-vector of clipping LIMITS [v_min,v_max,i_min,i_max] (see MDPWL)
    %   TOL - 1-3 vector of tolerance(s) {[rel],[tolv,toli],[v,i,rel]} (see MDPWL)

        % narginchk(0,5);

        % empty class constructor
        if nargin == 0; return; end
        if nargin < 2, ge = []; end
        if nargin < 3, tc = []; end

        if isa(source,'mdpwl') && size(source,1) == numel(ge) && size(source,2) == numel(tc)
        % construct from given set of IV curves
            obj = interpolantfromIVpps(source,ge,tc,varargin{:});
        elseif isstruct(source) && all(isfield(source,{'Iph_ref','Io_ref','Rsh_ref','Rs','G_ref'}))
        % construct from One Diode Model
            obj = interpolantfromODM(source,ge,tc,varargin{:});
        else
            error('ModuleInterpolant: Sorry, this only works with ODM so far'); 
        end
        

        function obj = interpolantfromODM(ODM,ge,tc,Lims,RelTol,varargin)
        % OBJ = INTERPOLANTFROMODM(ODM,GE,TC,LIMS,RELTOL)
        % See documentation in MODULEINTERPOLANT constructor

            % narginchk(1,5);
            
            % Load default evaluation grid, limits, and PWL approx. tolerances; if necessary
            if nargin < 2 || isempty(ge), ge = getSimOption('Gegrid'); end
            if nargin < 3 || isempty(tc), tc = getSimOption('Tcgrid'); end
            if nargin < 4, Lims = []; end
            if nargin < 5 || isempty(RelTol), RelTol = getSimOption('RelTol'); end

            % Using ge ~= 0 improves interpolation for low irradiances:
            % if any(ge == 0), ge = ge + getSimOption('minGHI'); end
            ge = sort(ge); tc = sort(tc);
                
            % check Module for required fields, check that all parameters are finite scalars
            ODM = checkODM(ODM);
            
            % Recommended limits [Vmin Vmax Imin Imax]
            RecLims(2) = onediodemodel2(translateODM(ODM,ge(end),tc(1)),0);    % max Voc
            RecLims(4) = onediodemodel(translateODM(ODM,ge(end),tc(end)),0);   % max Isc
            RecLims(3) = -RecLims(4);                                        % - Isc
            RecLims(1) = onediodemodel2(translateODM(ODM,0,25),RecLims(4));  % Vr @ max Isc, G = 0

            if isempty(Lims), Lims = RecLims; end
            assert(numel(Lims) == 4,'Limits do not seem right');

            if any([-1 1 -1 1].*(Lims - RecLims) < 0)
                RecLims([1,3]) = min(Lims([1,3]),RecLims([1,3]));
                RecLims([2,4]) = max(Lims([2,4]),RecLims([2,4]));
                printlim = @(x) sprintf('%0.1f to %0.1f V, %0.1f to %0.1f A',x);
                warning('Adjusting interpolant limits (%s) to recommended values (%s)',...
                    printlim(Lims),printlim(RecLims));
                Lims = RecLims;
            end

            [Tol,tolv,toli] = parsetolerance([RelTol*ODM.Vmp_ref,RelTol*ODM.Imp_ref,RelTol]);

            waitbar_obj = optwaitbar(0,'Preparing ODM Module Interpolant...');

            F(ge,tc);
            
            % Generate the set of reference IV curves by recursive refining of grids
            % (see functions F, G, E below)
            [x,IVcurveSet] = fitinterpolant(@F,{ge,tc},1,@G,@E);

            waitbar_obj.update(1.0,'Creating Interpolant from IV-curve set...');

            % Create IVppInterpolant object from reference curves
            obj = interpolantfromIVpps(IVcurveSet,x{:},ODM);
            obj.source = ODM;
            obj.info = struct('ge',ge,'tc',tc,'Lims',Lims,'RelTol',RelTol);
            obj.Tol = Tol;

            function V = F(g,t)
            % Value function for FITINTERPOLANT: return PWL curves for G, T
                n = numel(g);
                waitbar_obj.reset();
                for k = n:-1:1
                    pars = translateODM(ODM,g(k),t(k));
                    V(k) = onediodepwlapprox(pars,Tol,Lims,varargin{:});

                    prog = (n-k)/n;
                    waitbar_obj.update(prog,sprintf('Evaluating ODM at grid point %d/%d',n-k,n),'-addtime');
                end
                V = reshape(V,size(g));
            end

            function Q = G(g,t,V,Gq,Tq)
            % Interpolation function for FITINTERPOLANT
                pmi = ModuleInterpolant();
                pmi.Gref = g;
                pmi.Tref = t;
                pmi.IVpps = V;
                Q = getIVpp(pmi,Gq,Tq);
            end

            function e = E(Vf,Vg)
            % Error function for FITINTERPOLANT
                [Pmp0,Vmp0] = arrayfun(@mpp,Vf);
                [Pmp,~] = arrayfun(@mpp,Vg);
                tolp = min(tolv(Vmp0).*ODM.Imp_ref,toli(Pmp0./Vmp0).*ODM.Vmp_ref,'omitnan');
                e = abs(Pmp - Pmp0)./tolp;
            end
        end
        
        function obj = interpolantfromIVpps(IVpps,ge,tc,P)
        % OBJ = INTERPOLANTFROMIVPPS(IVPPS,GE,TC,NOMPARS)
        % See documentation in MODULEINTERPOLANT constructor

            narginchk(3,4);
            assert(numel(ge) > 1 && numel(tc) > 1,'GE and TC must have at least 2 values');
            assert(all(size(IVpps) == [numel(ge),numel(tc)]),'IVPPS does not match GE, TC');
            if nargin < 4 || isempty(P), P = struct(); end
            assert(isstruct(P),'Expecting nominal parameter structure PARAMS');
        
            obj = ModuleInterpolant(); % start with empty class
            obj.Gref = ge;
            obj.Tref = tc;
            obj.IVpps = IVpps;
            obj.source = 'IVpps';

            % set propperties to default values
            DEF.name = '';	
            DEF.G0 = 1000;
            DEF.Tc0 = 25;
            DEF.area = 0;
            DEF.geom = pvArea();
            
            % Put any non-recognized params into obj.info
            defnames = fieldnames(DEF)';
            obj.info = rmfield(P,intersect(fieldnames(P),defnames));
            P = rmfield(P,fieldnames(obj.info)); 
            P = completestruct(P,DEF); % complete with defaults
            for f = defnames, obj.(f{1}) = P.(f{1}); end

            % Create interpolants and get STC properties
            obj = obj.loadobj(obj);
            
            if obj.area == 0 && obj.geom.area > 0, obj.area = obj.geom.area; end
            if obj.area == 0
                warning('ModuleInterpolant:noarea',...
                    'Unknown area: MPPefficiency(obj...) will not work');
            end
        end
        
    end

    function IVpp = getIVpp(obj,g,t)
    % IVpp = GETIVPP(OBJ,G,T) - Get interpolated MDPWL curve for irradiance G and cell-temp T
    %   using bilinear interpolation according to Tsuno et al. 2009 [1]
    %
    %   G - effective (i.e. spectral/IAM-corrected, POA) irradiance (W/m²)
    %   T - cell-temperature (°C)
    %
    % [1] Y. Tsuno, Y. Hishikawa, and K. Kurokawa, “Modeling of the I–V curves of the PV modules
    % using linear interpolation/extrapolation,” Solar Energy Materials and Solar Cells,
    % vol. 93, no. 6–7, pp. 1070–1073, Jun. 2009.

        assert(isequal(size(g),size(t)),'Non-matching G,T');
        if ~isscalar(g)
            IVpp = arrayfun(@obj.getIVpp,g,t);
            return;
        end
        
        msg = 'beyond interpolation grid, returning curve for';
        if g > obj.Gref(end)
            msg = sprintf('Irradiance (%0.1f W/m²) %s %0.1f W/m², %0.1f°C',g,msg,obj.Gref(end),t);
            g = obj.Gref(end);
        elseif g < obj.Gref(1)
            msg = sprintf('Irradiance (%0.1f W/m²) %s %0.1f W/m², %0.1f°C',g,msg,obj.Gref(1),t);
            g = obj.Gref(1);
        end
        if t > obj.Tref(end)
            msg = sprintf('Temperature (%0.1f °C) %s %0.1f°C, %0.1f W/m²',t,msg,obj.Tref(end),g);
            t = obj.Tref(end);
        elseif t < obj.Tref(1)
            msg = sprintf('Temperature (%0.1f °C) %s %0.1f°C, %0.1f W/m²',t,msg,obj.Tref(1),g);
            t = obj.Tref(1);
        end
        if numel(msg) > 50
            warning('ModuleInterpolant:extrap',msg);
        end

        % find the indexes of the closest available IVpps
        [e1,e2] = getneighbours(g,obj.Gref);
        [t1,t2] = getneighbours(t,obj.Tref);

        % interpolate on temperature for each grid irradiance
        a = (t - obj.Tref(t1))/(obj.Tref(t2)-obj.Tref(t1));
        IVpp_e1tx = interpolate(obj.IVpps(e1,t1),obj.IVpps(e1,t2),a);
        IVpp_e2tx = interpolate(obj.IVpps(e2,t1),obj.IVpps(e2,t2),a);
        % ... and finally on irradiance
        b = (g - obj.Gref(e1))/(obj.Gref(e2)-obj.Gref(e1));
        IVpp = interpolate(IVpp_e1tx,IVpp_e2tx,b);
        IVpp.tol = obj.Tol;

        function [a,b] = getneighbours(x,P)
        % locate the indexes a and b of the closest and 2nd closest values in P to x
            a = find(x >= P,1,'last');
            if isempty(a), a = 1; b = 2;
            elseif a == numel(P), b = a-1; 
            elseif P(a+1) - x < x - P(a)
                b = a; a = a+1;
            else
                b = a+1;
            end
        end
    end

    function newobj = scale(obj,Sv,Si)
    % OBJ2 = SCALE(OBJ,Sv,Si) - Scale the IV-curve interpolant OBJ by Sv in V, and
    %   by Si in I. This can be used to represent arrays of obj elements (Sv,Si > 1 & integers),
    %   sub-elements of obj (1/Sv,1/Si > 1 & integers), and/or reverse connections (Sv,Si < 0).

        newobj = copy(obj); % start with original object
        newobj.name = sprintf('%s × %s scaled %s',nicescale(Sv),nicescale(Si),obj.name);
        newobj.area = obj.area/(Sv*Si);
        newobj.geom = pvArea(); % has no further meaning
        if isstruct(newobj.source)
            newobj.source.scale = [Sv,Si];
        elseif ischar(newobj.source)
            newobj.source = sprintf('%s × %s scaled %s',nicescale(Sv),nicescale(Si),newobj.source);
        else
            warning('Could not store Scale Factors on scaled ModuleInterpolant');
        end
        
        % Scale IVpps: v' = Sv·v, i' = Si·i
        newobj.IVpps = arrayfun(@(p) scale(p,Sv,Si),newobj.IVpps);

        % Recalculate interpolants and nominal values
        newobj = newobj.loadobj(newobj);
        
        function s = nicescale(x)
           if x >= 1, s = sprintf('%0.3g',x); else, s = sprintf('1/%0.3g',1/x); end
        end
    end

    function [Pmp,Vmp,Imp,c] = getMPP(obj,g,t,method)
    % [Pmp,Vmp,Imp,C] = getMPP(OBJ,G,T,METHOD) - return the Maximum-Power-Point of OBJ at
    %   irradiance G (W/m²) and cell-temperature T (°C).
    %
    %   G - scalar (constant) or array of effective irradiance values (W/m²).
    %     NOTE: Any G < 0 is treated as 0.
    %   T - scalar (constant) or array of cell-temperature values (°C)
    %
    %   METHOD - (optional) string. If METHOD is 'interp' or omitted, IV-curves are not generated
    %     or interpolated: pre-loaded Pmp(g,t) and Imp(g,t) interpolants are called directly.
    %     If method is 'PWL' or 'IVpps', the interpolated IV-curve is generated, and the MPP
    %     is searched over the true curve. 
    %
    % See also: MDPWL.MPP

        narginchk(3,4);
        if nargin < 4 || isempty(method), method = 'interp'; end

        [g,t] = compatiblesize(g,t);
        Pmp = zeros(size(g));
        Vmp = zeros(size(g));
        Imp = zeros(size(g));
        c = zeros(size(g));

        notdark = g > 0;
        if ~any(notdark), return; end

        switch lower(method)
            case {'mppinterp','interp','approx'}
                gx = single(g(notdark)); % double is required by interpolants
                tx = single(t(notdark));
                Pmp(notdark) = obj.PmpInterp(gx,tx); 
                Imp(notdark) = obj.ImpInterp(gx,tx);
                Vmp(notdark) = Pmp(notdark)./Imp(notdark);
                Vmp(Imp == 0) = 0;
                
                if nargin > 3
                   c(notdark) = obj.cInterp(gx,tx); 
                end
            case {'ivpp','ivpps','pwl','mdpwl'}
                
                DELTA_VMP = 0.01;
                x = 1+[-1,1]*DELTA_VMP; % [V-dV,V+dV]/V near Vmp, used to get C
                
                ok = g >= obj.Gref(1) & g <= obj.Gref(end) & t >= obj.Tref(1) & t <= obj.Tref(end);
                if any(~ok)
                   Pmp(~ok) = NaN;
                   Vmp(~ok) = NaN;
                   Imp(~ok) = NaN;
                   notdark = notdark & ok;
                end
        
                for j = find(notdark(:))'
                    pp = obj.getIVpp(g(j),t(j));
                    [Pmp(j),Vmp(j)] = mpp(pp);
                    Imp(j) = Pmp(j)/Vmp(j);
                    if nargout > 3
                        y = pp.val(Vmp(j)*x)/Imp(j);
                        c(j) = -(y(2)+y(1)-2)/DELTA_VMP^2; % y"
                    end
                end
            case {'odm','exact'}
                try 
                    % ODM = checkODM(obj.source);
                    ODM = obj.source;
                    translateODM(ODM,800,40);
                catch ERR
                    error('OBJ.source must be a working ODM structure: \n\n%s',getReport(ERR))
                end
                [~,Vmp] = getMPP(obj,g,t,'interp'); % use interpolated Vmp as seeds
                [Pmp(notdark),Vmp(notdark),Imp(notdark),~,~,c(notdark)] = ...
                    arrayfun(@(g,t,v0) onediodeMPP(translateODM(ODM,g,t),obj.Tol,v0),...
                    g(notdark),t(notdark),Vmp(notdark));
            otherwise
                error('ModuleInterpolant2:getMPP:method','Unknown method');
        end
    end

    function eta = MPPefficiency(obj,g,t,varargin)
    % E = APPROXMPPEFFICIENCY(OBJ,G,T,METHOD) - Get module MPP-efficiency at irradiance(s) G and
    %   cell-temperature(s) T. Optional METHOD is passed to GETMPP.
    
        eta = getMPP(obj,g,t,varargin{:})./(g*obj.area);
        eta = max(eta,0);
    end

    function Isc = getIsc(obj,g,t)
    % Get Isc at irradiance(s) g and cell-temperature(s) t, from scalar interpolant IscInterp
        dark = g <= 0;
        Isc = zeros(size(g));
        if all(dark), return; end
        Isc(~dark) = obj.IscInterp(double(g(~dark)),double(t(~dark)));
    end

    function Voc = getVoc(obj,g,t)
    % Get Voc at irradiance(s) g and cell-temperature(s) t, from scalar interpolant VocInterp
        dark = g <= 0;
        Voc = zeros(size(g));
        if all(dark), return; end
        Voc(~dark) = obj.VocInterp(double(g(~dark)),double(t(~dark)));
    end
    
    function varargout = InterpolationTest(obj,n)
    % OBJ.INTERPOLATIONTEST(N) - Evaluate the precision of MODULEINTERPOLANT object OBJ and plot
    %   the results. Namely, generate a mesh of irradiance and cell-temperature values that is
    %   N times finer than OBJ.Gref, OBJ.Tref, and compare Vmp(g,t), Pmp(g,t) accross this grid 
    %   when evaluated directly from the underlying ODM, from the PWL interpolation, and using 
    %   scalar interpolants.

        if nargin < 2 || isempty(n), n = getSimOption('meshref'); end

        gg = resamplestructure(obj.Gref(:),n);
        gg(1:n:end,:) = [];
        tt = resamplestructure(obj.Tref(:),n);
        tt(1:n:end,:) = [];

        [gg,tt] = meshgrid(gg,tt);
        pmp = zeros([3,size(gg)]);
        vmp = zeros([3,size(gg)]);

        w = stopwatch();
        fprintf('\nEvaluating Module Interpolant precision...\n\n');
        
        warning_resetter = naptime('checkODM:nomchk');  %#ok<NASGU>
        ODM = checkODM(obj.source);
        
        fprintf('\tCalculating from interpolants...');
        w.resetlocal();
        % Case 3: get MPP from direct interpolation structures
        [pmp(3,:),vmp(3,:)] = getMPP(obj,gg(:),tt(:),'interp');
        w.add2lap('interp');
        fprintf('\t(%0.4f seconds elapsed)\n',w.laps.interp);
        
        fprintf('\tCalculating from IVpps: 000000/000000');
        w.resetlocal();
        for j = 1:numel(gg)
        % Case 2: get complete IVpp curve, then search MPP
            [pmp(2,j),vmp(2,j)] = obj.getMPP(gg(j),tt(j),'ivpp');
            if mod(j,100) == 0
                fprintf([repmat('\b',1,13) '%06d/%06d'],j,numel(gg));
            end
        end
        w.add2lap('IVpp');
        fprintf('\t(%0.2f seconds elapsed)\n',w.laps.IVpp);
        
        fprintf('\tCalculating Pmp by ODM: 000000/000000');
        w.resetlocal();
        for j = 1:numel(gg)
        % Case 1: "Real" MPP (directly from ODM)
            % [pmp(1,j),vmp(1,j)] = obj.getMPP(gg(j),tt(j),'exact');  
            [pmp(1,j),vmp(1,j)] = onediodeMPP(translateODM(ODM,gg(j),tt(j)),obj.Tol,vmp(3,j));
            if mod(j,20) == 0
                fprintf([repmat('\b',1,13) '%06d/%06d'],j,numel(gg));
            end
        end
        w.add2lap('ODM');
        fprintf('\t(%0.2f seconds elapsed)\n\n',w.laps.ODM);

        errPmp = squeeze(pmp(2,:,:)-pmp(1,:,:))/obj.source.Pmp_ref;
        errVmp = squeeze(vmp(2,:,:)-vmp(1,:,:))/obj.source.Vmp_ref;
        errPmp2 = squeeze(pmp(3,:,:)-pmp(1,:,:))/obj.source.Pmp_ref;
        errVmp2 = squeeze(vmp(3,:,:)-vmp(1,:,:))/obj.source.Vmp_ref;

        errstats = @(x) sprintf('MBE = %0.2e, RMS = %0.2e, P95 = %0.2e/%0.2e',...
                                 mean(x(:)), std(x(:)), prctilew(x,2.5), prctilew(x,97.5));
        fprintf('\terr_Pmp/Pmp0: %s\n',errstats(errPmp));
        fprintf('\terr_Vmp/Vmp0: %s\n\n',errstats(errVmp));

        h = GUIfigure('ModuleInterpolant','Module Interpolant Test',[0.15 0.08 0.7 0.8]); clf(h);

        [g0,t0] = meshgrid(obj.Gref,obj.Tref);
        
        subplot(221); 
        contourf(gg,tt,errPmp,10,'linecolor','none'); colorbar();
        hold on; plot(g0(:),t0(:),'k+','MarkerSize',2);
        axesandtitle('| P_{mp}(IVpp) – P_{mp}(ODM) | / P_{mp}^{ 0}'); 

        subplot(222); 
        contourf(gg,tt,errPmp2,10,'linecolor','none'); colorbar();
        hold on; plot(g0(:),t0(:),'k+','MarkerSize',2);
        set(gca(),'Xscale','log');
        axesandtitle('| P_{mp}(getMPP) – P_{mp}(ODM) | / P_{mp}^{ 0}');

        subplot(223); 
        contourf(gg,tt,errVmp,10,'linecolor','none'); colorbar();
        hold on; plot(g0(:),t0(:),'k+','MarkerSize',2);
        axesandtitle('| V_{mp}(IVpp) – V_{mp}(ODM) | / V_{mp}^{ 0}'); 

        subplot(224); 
        contourf(gg,tt,errVmp2,10,'linecolor','none'); colorbar();
        hold on; plot(g0(:),t0(:),'k+','MarkerSize',2);
        set(gca(),'Xscale','log');
        axesandtitle('| V_{mp}(getMPP) – V_{mp}(ODM) | / V_{mp}^{ 0}');
        
        if nargout > 0, varargout = {pmp,vmp,w}; end

        function axesandtitle(titlestr)
            title(titlestr);
            xlabel('Ge (W/m²)'); 
            ylabel('Tc (°C)');
        end
    end
    
    function plot(ModIV,Tc,Ge,ax)
    % plot(OBJ,[Tc,Ge,ax]) - Plot module IV curves at cell-temperatures Tc and irradiance Ge
    % plot(OBJ,@(g,t,w) FCN,[Ge, ax]) - Plot module IV curves at cell-temperatures calculated
    %   from cell-temperaturemodel FCN, for irradiances Ge.
    %
    % See also: PLOTONEDIODE
    
        if nargin < 2, Tc = []; end
        if nargin < 3 || isempty(Ge)
            if numel(Tc) > 1, Ge = 1000; else, Ge = 200:200:1000; end
        end
        if isempty(Tc)
            if numel(Ge) > 5, Ge = 1000; else, [Ge,Tc] = meshgrid(Ge(:),[25,40]); end
        end
        if nargin < 4, ax = gca(); end
            
        if isnumeric(Tc)
            [Ge,Tc] = compatiblesize(Ge,Tc);
            title(ax,ModIV.name);
        elseif isa(Tc,'function_handle')
            Tc = Tc(Ge,25,1);
            title(ax,[ ModIV.name ' @ Ta = 25°C, WS = 1 m/s' ]);
        end
        set(get(ax,'title'),'interpreter','none');
        
        set(ax,'colororder',jet(numel(Ge))*0.9);
        arrayfun(@(g,t) plot(ModIV.getIVpp(g,t),ax),Ge,Tc);
        xlabel(ax,'Module voltage [V]'); 
        ylabel(ax,'Module current [A]');
        legend(ax,arrayfun(@(g,t) sprintf('%0.0f W/m^2, Tc %0.0f°C',g,t),Ge,Tc,'unif',0),...
            'location','southwest','EdgeColor',[1 1 1]*0.8);
        axis(ax,[0 ModIV.Voc0 0 ModIV.Isc0*1.2]);
    end
end

methods(Static)
    function obj = loadobj(obj)
    % obj = LOADOBJ(obj) - Generate STC parameters and interpolants on object load

        if ~isa(obj.IVpps,'mdpwl')
        % Update deprecated MODULEINTERPOLANT objects
            for j = numel(obj.IVpps):-1:1
                x = obj.IVpps{j}.breaks';
                n = numel(x);
                y = [obj.IVpps{j}.coefs(:,2);0];
                y(n) = obj.IVpps{j}.coefs(n-1,1)*(x(n)-x(n-1))+y(n-1);
                pp(j) = mdpwl(x,y,0);
            end
            obj.IVpps = reshape(pp,size(obj.IVpps));
        end
    
        if obj.G0 ~= 1000, warning('CAUTION: obj.G0 is not STC 1000 W/m^2 !!!'); end
        if obj.Tc0 ~= 25, warning('CAUTION: obj.Tc0 is not STC 25°C !!!'); end

        % Get STC model parameters
        [~,obj.Vmp0,obj.Imp0] = getMPP(obj,obj.G0,obj.Tc0,'pwl');
        pp = obj.getIVpp(obj.G0,obj.Tc0);
        obj.Isc0 = pp.val(0);
        obj.Voc0 = pp.inval(0);

        % Use a finer mesh (1:1/MESHREF:N) for scalar interpolants
        MESHREF = getSimOption('meshref');
        ge = interp1(obj.Gref,1:1/MESHREF:numel(obj.Gref));
        tc = interp1(obj.Tref,1:1/MESHREF:numel(obj.Tref));

        if isempty(obj.cInterp)
        % Evaluate MPP, Voc, and Isc of each curve
        
            if str2double(regexprep(version(),'(\d+\.\d+).*','$1')) <= 9.2
                METHOD = 'pchip';
            else
                METHOD = 'makima';
            end
        
            % [Pmp,Vmp,Imp,c] = getMPP(obj,ge',tc,'ivpp');
            
            DELTA_VMP = 0.01;
            x = 1+[-1,1]*DELTA_VMP; % [V-dV,V+dV]/V near Vmp, used to get C
        
            for i = numel(ge):-1:1
                for j = numel(tc):-1:1
                    pp = obj.getIVpp(ge(i),tc(j));
                    [Pmp(i,j),Vmp(i,j)] = pp.mpp();
                    Isc(i,j) = pp.val(0);
                    Voc(i,j) = pp.inval(0);
                    if Pmp(i,j) > 0
                        Imp(i,j) = Pmp(i,j)./Vmp(i,j);
                        y = pp.val(Vmp(i,j)*x)/Imp(i,j);
                        c(i,j) = -(y(2)+y(1)-2)/DELTA_VMP^2; % y"
                    % else Imp(i,j) = c(i,j) = 0
                    end
                end
            end
            
            if isstruct(obj.source) && all(isfield(obj.source,{'nVth0','Iph0','muIsc','Io0','Rsh0','Rs','G0'}))
            % Refine MPP values using original ODM
                [ge,tc] = ndgrid(ge,tc);
                [Pmp,Vmp,Imp,Voc,Isc,c] = arrayfun(@(g,t,v) onediodeMPP(translateODM(obj.source,g,t),obj.Tol,v),ge,tc,Vmp);
                ge = ge(:,1); tc = tc(1,:);
            end
            
            % Create interpolants for each variable
            obj.PmpInterp = griddedInterpolant({ge,tc},single(Pmp),METHOD);
            obj.VocInterp = griddedInterpolant({ge,tc},single(Voc),METHOD);
            obj.IscInterp = griddedInterpolant({ge,tc},single(Isc),METHOD);
            obj.cInterp = griddedInterpolant({ge,tc},single(c),METHOD);
            
            % Imp yields a much smoother interpolant than Vmp
            %Imp = Pmp./Vmp; 
            Imp(Vmp == 0) = 0;
            obj.ImpInterp = griddedInterpolant({ge,tc},single(Imp),METHOD);
        end
    end
end

end
