	classdef arraysolverplotter < plotter
	properties
        subplothandles  % structure of sub-plot handles {data,shades,ivpp,pvpp}
		axislimits      % structure with axis limits [xmin xmax ymin ymax] for each sub-plot
    end 
    properties (Hidden = true)
		trackeroutline  % Polygon object with tracker outline
		inverterwindow  % Polygon object for inverter window
        textcoords      % last written text box in subplot(1,1)
    end
    methods
        function obj = arraysolverplotter(plotting,trckoutline,invlimits)
        % Draw one of the two instances of the arraysolverplotter:
        %   - Waitbar, with option button to switch to Plotter
        %     obj = arraysolverplotter() or arraysolverplotter(false,...)
        %   - Plotter: 4 axes plot with control buttons, and option to switch to waitbar  
        %     obj = arraysolverplotter(true,trckoutline,invlimits)
        %           trckoutline: tracker outline polygon object
        %           invlimits: [Vmin Vmax Ps0 Pmax]
        
            UI = runningfromUI();
            if nargin < 1, plotting = UI; end

            % Dummy inputs for testing
            if nargin < 2, trckoutline = polygon(1,1); end
            if nargin < 3, invlimits = [0 Inf 0 Inf]; end
            
            % Initialize base plotter object (it will call this class's drawplotter/drawwaitbar)
            obj@plotter(plotting,...
                {0,'Initializing...','Name','Electrical Solver'},...
                {trckoutline,invlimits});
 
            if ~obj.plotting
            % Try to set sub-plot axis-limits, obj.trckoutline and obj.inverterwindow
            % ...also when ~plotting, never know when switchtoplotter will be called
                setaxislimits(obj,trckoutline,invlimits);
            end
        end

        function drawplotter(obj,trckoutline,invlimits)
            
            if isempty(obj.axislimits)
                assert(nargin == 3,'First call to DRAWPLOTTER requires 3 arguments');
                setaxislimits(obj,trckoutline,invlimits);
            end

            drawplotter@plotter(obj,'Shading-Analysis Plotter','big');

            obj.subplothandles = struct();
			obj.subplothandles.data = subplot(221);
			obj.subplothandles.shades = subplot(222);
			obj.subplothandles.ivpp = subplot(223);
			obj.subplothandles.pvpp = subplot(224);
			
			% Initialize blackboard
            set(obj.subplothandles.data,'Visible','off')
            obj.textcoords = [0,1];
                        
			% Initialize beam-shading subplot
			subplot(obj.subplothandles.shades)
			xlabel('Tracker x');
			ylabel('Tracker y');
			set(gca,'XTickLabel',[]);
			set(gca,'YTickLabel',[]);
			axis(obj.axislimits.shades)
			axis square
            set(gca,'visible',false);

			% Initialize array-IV-curve subplot
			subplot(obj.subplothandles.ivpp)
			xlabel('Array Voltage [V]')
			ylabel('Array Current [A]')
			axis(obj.axislimits.ivpp)
            set(gca,'yscale','log');
            % xlim(obj.axislimits.ivpp(1:2));
			axis square

			% Initialize array-PV-curve subplot
			subplot(obj.subplothandles.pvpp)
			xlabel('Array Voltage [V]')
			ylabel('Array Power [kW]')
			axis(obj.axislimits.pvpp) 
            set(gca,'yscale','log');
            % xlim(obj.axislimits.ivpp(1:2));
			axis square;

			reset(obj);
        end
                
        function setaxislimits(obj,trckoutline,Inverter)
        % Set sub-plot axis limits, based on trckoutline and invlimits, if available
        % set also obj.trckoutline and obj.inverterwindow properties
        %   trckoutline: tracker border polygon (polygon object)
        %   invlimits: inverter operational envelope [Vmin Vmax Ps0 Pmax]
        
            obj.axislimits.data = [0 1 0 1];

            if nargin > 1 && isa(trckoutline,'polygon') && trckoutline.area > 0
                obj.trackeroutline = trckoutline;
                [w,h] = rectangleproperties(trckoutline);
                obj.axislimits.shades = [-0.5 0.5 -0.5 0.5]*1.2*max(w,h);
            else
                obj.trackeroutline = polygon(); 
            end
            
            if nargin > 2 && isstruct(Inverter) 
                Plim = [Inverter.Ps0,Inverter.Pdc0];
                Vlim = [Inverter.MPPTLow,Inverter.MPPTHi];
                Ilim = [0,Inverter.IdcMax];

                % Inverter window (on IV space) is the intersection of P/V limits and I/V limits:
                vv = linspace(Vlim(1),Vlim(2),20);
                p = polygon([vv,fliplr(vv)]',[Plim(2)./vv,Plim(1)./fliplr(vv)]');
                p = intersectpolygons(p,polygon(Vlim,Ilim));
                p = refinepoly(p,1/40,@(x,y) abs(x/diff(Vlim)));

                obj.inverterwindow = p;
                obj.axislimits.ivpp = [0 Vlim(2)*1.1 1e-1 Ilim(2)];
                obj.axislimits.pvpp = [0 Vlim(2)*1.1 1e-1 Plim(2)/1000*1.1]; 
            else
                obj.inverterwindow = polygon();
                obj.axislimits.ivpp = [0 Inf 1e-1 Inf];
                obj.axislimits.pvpp = [0 Inf 1e-1 Inf];
            end
        end
		
		function reset(obj,axeshandles)
        % reset(obj), or reset(obj,[hi,hj...])
        % clear the contents of all/selected sub-plot axeshandles
        
            assert(ishandle(obj.fighandle),'arraysolverplotter:noobj','Closed by user');
			%figure(obj.fighandle);
			if nargin < 2 || isequal(axeshandles,':')
				axeshandles(1) = obj.subplothandles.data;
				axeshandles(2) = obj.subplothandles.shades;
				axeshandles(3) = obj.subplothandles.ivpp;
				axeshandles(4) = obj.subplothandles.pvpp;
			end
			
			% Clear blackboard (data subplot)
			if any(axeshandles == obj.subplothandles.data)
				subplot(obj.subplothandles.data)
				hold on; cla
                obj.textcoords = [0,1];
			end

			% Clear beam-shading subplot
			if any(axeshandles == obj.subplothandles.shades)
				subplot(obj.subplothandles.shades)
				hold on; cla
% 				polyplot(obj.trackeroutline,[0.2 0.2 0.6]);
			end
			
			% Clear array-IV-curve subplot
			if any(axeshandles == obj.subplothandles.ivpp)
				subplot(obj.subplothandles.ivpp)
				hold on; cla
                polyplot(obj.inverterwindow,'none',1,'linestyle',':');
			end
			
			% Clear array-PV-curve subplot
			if any(axeshandles == obj.subplothandles.pvpp)
				subplot(obj.subplothandles.pvpp)
				hold on; cla
                p = obj.inverterwindow;
                p.y = p.x.*p.y/1000;
				polyplot(p,'none',1,'linestyle',':');
			end
		end
		
		function plotarraycurves(obj,IVpp,varargin)
			
            assert(ishandle(obj.fighandle),'arraysolverplotter:noobj','Closed by user');
			if nargin < 3 || isempty(varargin), varargin = {'b-'}; end
            
            pp = copy(IVpp);
            pp.clip(obj.axislimits.ivpp);

			figure(obj.fighandle)
			subplot(obj.subplothandles.ivpp)
			plot(pp.x,pp.y,varargin{:});

			subplot(obj.subplothandles.pvpp)
            xx = [pp.x];
            xx = unique([xx;linspace(0,max(xx),100)'],'sorted');
			plot(xx,xx.*pp.val(xx)/1000,varargin{:});
		end
		
		function plotpointoncurves(obj,Vx,Px,varargin)
			
            assert(ishandle(obj.fighandle),'arraysolverplotter:noobj','Closed by user');
			if nargin < 3 || isempty(varargin), varargin = {'ro'}; end
			
			figure(obj.fighandle)
			subplot(obj.subplothandles.ivpp)
			plot(Vx,Px./Vx,varargin{:});
		
			subplot(obj.subplothandles.pvpp)
			plot(Vx,Px/1000,varargin{:});
        end
		
        function addtextline(obj,txt)
            assert(ishandle(obj.fighandle),'arraysolverplotter:noobj','Closed by user');
            subplot(obj.subplothandles.data)
            txtid = text(obj.textcoords(1),obj.textcoords(2),txt);
            h = get(txtid,'Extent');
            obj.textcoords = obj.textcoords - [0,h(4)];
        end
        
		function plotshadedtrackers(obj,shadepoly,wt)
			
            assert(ishandle(obj.fighandle),'arraysolverplotter:noobj','Closed by user');
			nt = numel(shadepoly);
			subplot(obj.subplothandles.shades)
			[w,h,~] = rectangleproperties(obj.trackeroutline);
            nc = ceil(sqrt(nt));
            nr = ceil(nt/nc);
			dx = (0.5:nc)*1.2*w;
            dy = (0.5:nr)*1.2*h;
            [dx,dy] = meshgrid(dx,dy);
			obj.axislimits.shades = [0 nc*w 0 nr*h]*1.2;
			axis(obj.axislimits.shades); axis square;
			
            DARK_COLOR = [0.1 0.1 0.5]; % (background) dark blue
			LIGHT_COLOR = @(w) [0.2 0.2 0.8 w]; % alpha-weighted light blue
						
			hold on; cla
			for j = 1:nt
                polyplot(movepoly(obj.trackeroutline,dx(j),dy(j)),DARK_COLOR); % background
				if isnumeric(shadepoly{j})
					if shadepoly{j} == -1
                    % fully shaded, plot nothing 
					else
						polyplot(movepoly(obj.trackeroutline,dx(j),dy(j)),LIGHT_COLOR(1));
					end
                else
                    for i = 1:numel(shadepoly{j})
                        p = polygon.unpack(shadepoly{j}{i});
                        polyplot(movepoly(p,dx(j),dy(j)),LIGHT_COLOR(wt{j}(i)),DARK_COLOR);
                    end
				end
			end
			function q = movepoly(p,x,y)
				q = p;
				for k = 1:numel(q)
					q(k).x = q(k).x+x;
					q(k).y = q(k).y+y;
				end
			end
        end
        
        function exportfigure(obj,filename,resolution)
            
            assert(ishandle(obj.fighandle),'arraysolverplotter:noobj','Closed by user');
            if nargin < 3, resolution = 300; end
            figure(obj.fighandle)
            if isempty(dir('./fig')), mkdir('fig'); end
            print('-dpng',sprintf('-r%d',resolution),filename)
% 			export_fig(sprintf('./Frames/scen%s',filename),specs{:},obj.subplothandles.shades);
% 			export_fig(sprintf('./Frames/proj%s',filename),specs{:},obj.subplothandles.ivpp);
% 			export_fig(sprintf('./Frames/poly%s',filename),specs{:},obj.subplothandles.pvpp);
            cd('../')
        end
	end
end
