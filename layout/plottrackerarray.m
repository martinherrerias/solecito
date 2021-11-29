function varargout = plottrackerarray(Trck,varargin)
% PLOTTRACKERARRAY(TRCK) - Plot a representation of PV-plant layout TRCK.
% PLOTTRACKERARRAY(TRCK,az,el) - use sun-position(s) az,el. Defaults are the sunpaths at solstices.
% PLOTTRACKERARRAY(..,AX) - plot in axes handle AX
%
% PLOTTRACKERARRAY(..,'labels',X) - for X in {'on','off','all', NMAX} controls if (or how many)
%   mount ID labels are included in the plot. Default is 'on', which sets NMAX = 200.
%
% PLOTTRACKERARRAY(..,'idx',N,['axes',true]) - plot rotated coordinate axes for mount N. 
%   Default 'idx' chooses a mount near the center of the layout.
%
% PLOTTRACKERARRAY(..,'idx',N,..,'mask',true) - color mounts based on the shading mask for mount N,
%   Default is 'mask',true if an explicit idx is provided.
%
% PLOTTRACKERARRAY(..,'view',[AZ,EL]) - set view([AZ,EL])
%
% TRCK (Required fields):
%     TRCK.type - {1aF,1aC,1aV,2a,0a} and, depending on TRCK.type: { tilt, azimuth, slope, ..
%       tracklimits, backtracking, and groundcoverratio }. See MOUNTROTATIONS for details.
%     TRCK.origin - [longitude °N,latitude °E, elevation] of project origin
%     TRCK.rotation - Degrees N2W (CCW from North) of project's Y axis.
%     TRCK.geom.border - polygon object with tracker outline in Tracker-Coordinate-System.
%     TRCK.centers - 3·Ntr array of center point coordinates, for all trackers Ntr in the array.
%     TRCK.centerheight - scalar, required for representation of ground mesh
%   
%   Optional:
%     TRCK.tidx - Ntr·1 vector of integer-labels for each mount, default is (1:Ntr)'
%     TRCK.masks - Nat·Ntr or Ntr·Ntr boolean array
%     TRCK.axisoffset - vector from center of rotation to center of tracker, in TCS ([0;0;0])
%     TRCK.analysedtrackers - vector of Nat indices for trackers under analysis. (1:Ntr if omitted)
%
% See also: SAMPLETRACKERS, APPROXGROUNDMESH, MOUNTROTATIONS

    COLORS.sun = [1 0.8 0];
    COLORS.active = [0.3,0.3,0.8,1];
    COLORS.passive = [0.6,0.6,0.6,0.1];
    COLORS.gnd = [0,0.5,0.4;1,1,0.4];
    
    if ~isempty(varargin) && isscalar(varargin{end}) && ...
            isa(varargin{end},'matlab.graphics.axis.Axes')
        ax = varargin{end};
        varargin(end) = [];
    else
        ax = [];
    end

    opt.az = [];
    opt.el = [];
    opt.gnd = [];
    opt.shading = true;
    opt.azconv = 'N2E';
    opt.labels = 'off';
    opt.idx = [];
    opt.mask = 0;
    opt.axes = 1;
    opt.view = [80,20];
    opt.btn = [];
    opt = getpairedoptions(varargin,opt,'dealrest',4);
    
    if isempty(opt.az) || isempty(opt.el)
        [az,el] = sunpaths(Trck.origin);
    else
        [az,el] = deal(opt.az,opt.el);
        validateattributes(az,{'numeric'},{'real','vector','nonempty'});
        validateattributes(el,{'numeric'},{'real','vector','nonempty','size',size(az)});
        az = solarposition.fixazimuth(az,opt.azconv,'N2E');
    end
        
    parsestruct(Trck,{'centers'},'numeric','real','size',[3,NaN]);
    parsestruct(Trck,{'geom'},'class','pvArea','scalar');
    parsestruct(Trck,{'origin'},'numeric','real','size',[1,3]);
    
    Ntr = size(Trck.centers,2); 
    if isfield(Trck,'analysedtrackers')
        analysed = Trck.analysedtrackers;
        Nu = numel(Trck.analysedtrackers);
	else
		analysed = 1:Ntr;
        Nu = Ntr;
    end

    if isempty(opt.mask), opt.mask = ~isempty(opt.idx); end
    if ~isfield(Trck,'masks'), opt.mask = 0; end
    validateattributes(opt.mask,{'numeric','logical'},{'scalar','binary'});
    validateattributes(opt.axes,{'numeric','logical'},{'scalar','binary'});
    
    switch opt.labels
        case 'on', opt.labels = 200;
        case 'off', opt.labels = 0;
        case 'all', opt.labels = Inf;
        otherwise
            validateattributes(opt.labels,{'numeric'},{'integer','<=',Ntr});
    end

    % Draw trackers at start positions, store handles
    [F,updatefcn,sunvec,R] = plantlayout(Trck,az,el,'N2E');
    
    if ~isempty(opt.gnd) && (isa(opt.shading,'function_handle') || isempty(opt.shading))
        [Gnd,shfcn] = deal(opt.gnd,opt.shading);
    elseif opt.shading
        [shfcn,Gnd] = groundshading(Trck,sunvec,updatefcn);
    else
        Gnd = approxgroundmesh(Trck);
        shfcn = [];
    end
    
    if isnumeric(updatefcn)
        V0 = updatefcn;
    else
        V0 = updatefcn(1);
    end
        
    [scale,origin] = MinBoundCircle(Gnd.Points(:,1:2));
    
    if isempty(opt.idx)
        [~,opt.idx] = pdist2(Trck.centers(1:2,analysed)',origin,'euclidean','Smallest',1);
    else
        validateattributes(opt.idx,{'numeric'},{'scalar','integer','positive','<=',Nu});
    end
    origin = Trck.centers(:,analysed(opt.idx));
 
    if ~isempty(ax), axes(ax);
    else
        hfig = GUIfigure('trckfield'); 
        clf(hfig);
        set(hfig,'Name','Mount Layout'); set(hfig,'Numbertitle','off');
        ax = gca();
    end
    cla(ax);
    centeraxis(ax);
    ax.Color = 'none';
	hold(ax,'on');
    
    % Colors
    C = repmat(COLORS.passive,Ntr,1);
    C(analysed,:) = repmat(COLORS.active,numel(analysed),1);
    	
    using_masks = opt.mask;
    if using_masks
    % flag masked mounts for j'th analysed mount with green
        C(analysed(opt.idx),1:3) = [1,0,0];
        C(Trck.masks(opt.idx,1:3),2) = 1;
    end
    
    htr = patch(ax,'Faces',F,'Vertices',V0,'FaceVertexCData',C(:,1:3),'facecolor','flat',...
        'edgecolor','w','edgealpha',0.1,'FaceVertexAlphaData',C(:,4));

    % Plot ground mesh
	% colormap('summer');
    % trisurf(Gnd,'facecolor','interp','edgecolor','w','edgealpha',0.1,'SpecularStrength',0);
    
    clim = [min(Gnd.Points(:,3))-1,max(Gnd.Points(:,3))+1];
    z = mean(reshape(Gnd.Points(Gnd.ConnectivityList,3),[],3),2);
    GndColor = interpmatrix(z,linspace(clim(1),clim(2),size(COLORS.gnd,1))')*COLORS.gnd;
	hgnd = patch(ax,'Faces',Gnd.ConnectivityList,'Vertices',Gnd.Points,...
        'FaceVertexCData',GndColor,'facecolor','flat','edgecolor','none'); % 'w','edgealpha',0.1);

	% plot tracker centers and labels
    if opt.labels > 0
        if Ntr < opt.labels, s = 1:Ntr;  else, s = randperm(Ntr,opt.labels); end
        plot3(ax,Trck.centers(1,s),Trck.centers(2,s),Trck.centers(3,s),'r.')
        if ~isfield(Trck,'tidx'),Trck.tidx = (1:Ntr)'; end
        n0 = mean(R(:,3,:,:),3:4);
        s0 = Trck.centers + (Trck.axisoffset(3) + 0.5)*n0/norm(n0);
        text(ax,s0(1,s)',s0(2,s)',s0(3,s)',num2str(Trck.tidx(s)),'fontsize',8);
    end

    % % Print
    % view(180,90)
    % rez = 600; 
    % figpos = getpixelposition(hfig);
    % resolution = get(0,'ScreenPixelsPerInch');
    % set(hfig,'paperunits','inches','papersize',figpos(3:4)/resolution);...
    % set(hfig,'paperposition',[0 0 figpos(3:4)/resolution]);
    % print(hfig,'layout.png','-dpng',['-r',num2str(rez)],'-opengl')
    
    if opt.axes
        r0 = max(pdist(poly2vef(Trck.geom.border)));
        mountaxes = [0,0,0;1,0,0;nan(1,3);0,1,0;0,0,0;0,0,1]'*max(scale/2,r0);
        if isfield(Trck,'axisoffset'), mountaxes = mountaxes + Trck.axisoffset(:); end
        hax = plot3(ax,mountaxes(1,:),mountaxes(2,:),mountaxes(3,:),'r-');
    else
        hax = []; 
    end
    
    [spaz,spel] = sunpaths(Trck.origin); 
    sunpath = sph2cartV(90 - spaz - Trck.rotation,spel);

    % Draw sun-path and (some) sun vector
    hsunpath = plot3(ax, sunpath(:,1)*scale + origin(1),...
                         sunpath(:,2)*scale + origin(2),...
                         sunpath(:,3)*scale + origin(3),'-','color',COLORS.sun);
          
	hsun = plot3(ax,[origin(1),sunvec(1,1)*scale + origin(1)],...
                 [origin(2),sunvec(1,2)*scale + origin(2)],...
                 [origin(3),sunvec(1,3)*scale + origin(3)],...
                 '-','color',COLORS.sun); %draw vector
             
    % Rotate trackers to position (az,el) with max(el)
    last_t = 0;
    last_idx = 0;
    rotatetrackers(find(max(el)==el,1,'first'),opt.idx);
    Nt = numel(el);

    xlabel(ax,'x'); 
    ylabel(ax,'y');
    view(ax,opt.view);
	axis(ax,'equal');
    % set(hfig, 'toolbar','figure');
    ax.Clipping = 'off';
    ax.ZLimMode = 'manual';
    hold(ax,'off');

    if isempty(opt.btn)
        btn = uicontrol('Parent',ax.Parent,'Style','pushbutton','Visible','on',...
          'Units','normalized','Position',[1 0 0.8 0;0 1 0 0;0 0 0.2 0;0 0 0 0.1]*ax.Position');
    elseif ishandle(opt.btn) && isa(opt.btn,'matlab.ui.control.UIControl')
        btn = opt.btn;
    else
        btn = []; 
    end
    if ~isempty(btn)
        btn.String = 'Animate';
        btn.Callback = @animate;
        btn.UserData.running = false; % external kill switch
        if Nt > 1
            btn.Enable = 'on';
        else
            opt.btn.Enable = 'off';
        end
    end
    
    if nargout > 0, [varargout{1:2}] = deal(ax,@rotatetrackers); end
    
    function rotatetrackers(t,newidx,offset,masks,gndshade,trib)
        
        if nargin > 3 && ~isempty(masks), using_masks = masks; end
        if nargin < 3 || isempty(offset), offset = [0;0;0]; end
        
        if nargin < 2 || isempty(newidx), newidx = last_idx; end
        if last_idx ~= newidx
            last_idx = newidx;
            origin = Trck.centers(:,analysed(last_idx));
            
            sun = sunpath'*scale + origin;
            [hsunpath.XData,hsunpath.YData,hsunpath.ZData] = ...
                deal(sun(1,:),sun(2,:),sun(3,:));
            
            if using_masks
                C = repmat(COLORS.passive(1:3),Ntr,1);
                C(analysed,:) = repmat(COLORS.active(1:3),numel(analysed),1);
                C(analysed(last_idx),:) = [1,0,0];
                C(Trck.masks(opt.idx,:),2) = 1;
                htr.FaceVertexCData = C;
            end
        end

        if t ~= last_t
            if ~isnumeric(updatefcn)
                V = updatefcn(t);
                set(htr,'Vertices',V);
            end

            if ~isempty(shfcn)
                if nargin < 6
                    [gndshade,trib] = shfcn(t);
                end
                set(hgnd,'FaceVertexCData',GndColor.*(0.2+0.8*(1-gndshade).*trib));
            end

            if ishandle(hax)
                if size(R,4) > 1, tidx = t; else, tidx = 1; end
                if size(R,3) > 1, midx = analysed(last_idx); else, midx = 1; end
                pts = R(:,:,midx,tidx)*(mountaxes+offset) + origin;
                [hax.XData,hax.YData,hax.ZData] = deal(pts(1,:),pts(2,:),pts(3,:)); 
            end
            
            sun = sunvec(t,:)'*scale + origin;
            [hsun.XData(2),hsun.YData(2),hsun.ZData(2)] = deal(sun(1),sun(2),sun(3));
        end

        last_t = t;
    end
    
    function animate(src,~)
        
        if src.UserData.running
            src.String = 'Animate';
            src.UserData.running = false;
            return;
        else
            src.String = 'Stop';
            src.UserData.running = true;
        end
        
        for j = last_t + (1:Nt)
            t = mod(j-1,Nt)+1;
            if ~src.UserData.running, break; end
            tic();
            rotatetrackers(t);
            drawnow();
            dt = toc();
            pause(max(1/20-dt,0));
        end
        src.String = 'Animate';
        src.UserData.running = false;
    end
end

function [az,el] = sunpaths(origin)
% Get sunpaths (sunrise to sunset), at solstices

    lat = origin(2);

    STEP = 2;
    DECL = [-23.45,23.45];
    
    w0 = solarposition.sunsetangle(lat, DECL,'-finite');
    DECL(w0 == 0) = []; 
    w0(w0 == 0) = []; 
    m = ceil(2*w0/STEP);
    w = arrayfun(@(w,m) [linspace(-w,w,m),NaN],w0,m,'unif',0);
    w = cat(2,w{:});
    d = repelem(DECL,m+1);

    [az,el] = solarposition.eq2hor(lat,d,w,'N2E');
    az = az(:); el = el(:);
end


