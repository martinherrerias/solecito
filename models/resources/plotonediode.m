function varargout = plotonediode(pars,varargin)
% h = PLOTONEDIODE(PAR) - plot a PWL representation of one-diode-model 5-8 parameter structure
%   (or array of structures) PAR, as returned by TRANSLATEODM.
% h = PLOTONEDIODE(ODM) - equivalent to PAR = translateODM(ODM,Ge,Tc) for an automatic grid of
%   irradiance and cell-temperature values Ge, Tc.
%
% h = PLOTONEDIODE(..,'Lim',LIM) - pass limits (2/4-vector) to ONEDIODEPWLAPPROX, default [0,Voc]
% h = PLOTONEDIODE(..,'Tol',TOL) - pass tolerance (2/4-vector) to ONEDIODEPWLAPPROX
% h = PLOTONEDIODE(ODM,..,'Ge',G) - set custom grid of irradiance values
% h = PLOTONEDIODE(ODM,..,'Tc',T) - set custom grid of temperature values (compatible with Ge)
% h = PLOTONEDIODE(ODM,..,'Tc',@(g,t,w) FCN) - plot ODM curves at cell-temperatures calculated
%                                              from cell-temperaturemodel FCN, for irradiances Ge.
% h = PLOTONEDIODE(..,AX) - plot on axis AX

    [opt,varargin] = getpairedoptions(varargin,{'Lim','Tol','Ge','Tc'},{[],1e-3,[],[]});
    if ~isempty(varargin) && ishandle(varargin{end}) && isa(varargin{end},'matlab.graphics.axis.Axes')
        ax = varargin{end}; 
        varargin(end) = [];
    else, ax = gca();
    end

    parsestruct(pars,{},'opt',{'Iph','Io','Rs','Rsh','nVth',...
        'nDiode','Iph_ref','muIsc','Io_ref','Rsh_ref','Rs','Ns','area'},'-n','-r','-c');
     
    if all(isfield(pars,{'Iph','Io','Rs','Rsh','nVth'}))
        h = odmplot(pars);
    else
        try
            ODM = checkODM(pars);
        catch
            error('Expecting translated ODM parameters or full ODM structure');
        end
        h = plotfullODM(ODM,opt.Tc,opt.Ge);
    end
    if nargout > 0, varargout{1} = h; end
        
    function h = odmplot(P)
        set(ax,'colororder',jet(numel(P))*0.9);
        h = plot(arrayfun(@(p) onediodepwlapprox(p,opt.Tol,opt.Lim),P),varargin{:},ax);
        xlabel(ax,'voltage [V]'); 
        ylabel(ax,'current [A]');
        % axis(ax,[0 Inf 0 Inf]);
    end

    function h = plotfullODM(ODM,Tc,Ge)
    % plot(OBJ,[Tc,Ge,ax]) - Plot module IV curves at cell-temperatures Tc and irradiance Ge
    % plot(OBJ,@(g,t,w) FCN,[Ge, ax]) - Plot module IV curves at cell-temperatures calculated
    %   from cell-temperaturemodel FCN, for irradiances Ge.
    %
    % See also: PLOTONEDIODE

        if isempty(Ge)
            if numel(Tc) > 1, Ge = 1000; else, Ge = 200:200:1000; end
        end
        if isempty(Tc)
            if numel(Ge) > 5, Ge = 1000; else, [Ge,Tc] = meshgrid(Ge(:),[25,40]); end
        end
        if ~isfield(ODM,'name'), ODM.name = 'ODM'; end
        if isnumeric(Tc)
            [Ge,Tc] = compatiblesize(Ge,Tc);
            title(ax,ODM.name);
        elseif isa(Tc,'function_handle')
            Tc = Tc(Ge,25,1);
            title(ax,[ ODM.name ' @ Ta = 25°C, WS = 1 m/s' ]);
        end
        set(get(ax,'title'),'interpreter','none');

        P = translateODM(ODM,Ge,Tc);
        h = odmplot(P);

        legend(ax,arrayfun(@(g,t) sprintf('%0.0f W/m^2, Tc %0.0f°C',g,t),Ge,Tc,'unif',0),...
            'location','southwest','EdgeColor',[1 1 1]*0.8);
    end
end