classdef shadeanalysisplotter < plotter
    properties
        diffuseshading  % true: include POV-projections (wire and patches)
        groundshading
        
        axislimits      % [xmin xmax ymin ymax zmin zmax] for layout scene plotting
		axisoffset      % [x0,y0,z0] for plotting of moving axes and sun-vector
        subplothandles  % structure of sub-plot handles {scene,beamsh,povproj,shregions}
        
        Mounts
        sunaz
        sunel
        Gnd
    end 
    properties (Hidden = true)   
        polyh           % structure of object handles {trckaxis,sunvec,pjtrackers,pjshades}
        colors
        updatefcn = [];  % function handle, as returned by PLOTTRACKERARRAY
    end
    methods
        function obj = shadeanalysisplotter(plotting,varargin)
        % obj = shadeanalysisplotter(false)
        % obj = shadeanalysisplotter(plotting,Mounts,sunaz,sunel,Gnd,options)
        % 
        % Draw one of the two instances of the shadeanalysisplotter:
        %   - Waitbar, with option button to switch to Plotter
        %   - Plotter: 2/4 axes plot with control buttons, and option to switch to waitbar
        
            narginchk(6,6);
            
            % Initialize base plotter object (it will call this class's drawplotter/drawwaitbar)
            obj@plotter(plotting,'Shading Analysis',varargin);
            
            if nargin < 6 || any(cellfun(@isempty,varargin))
               assert(~plotting,'Plotting-enabled object requires full constructor call');
               setbtn(obj,[],'Enable','off')
            end
            varargin(end+1:5) = {[]};
            [obj.Mounts,obj.sunaz,obj.sunel,obj.Gnd,opt] = deal(varargin{:});
            if ~isempty(opt)
                obj.diffuseshading = opt.diffuseshading;
                obj.groundshading = opt.groundshading;
            end

            obj.axislimits = zeros(1,6);
            obj.axisoffset = [0;0;0];            
            obj.polyh = struct('trckaxis',NaN,'sunvec',NaN,'pjtrackers',NaN,'pjshades',NaN);
            obj.colors = getcolors();
            
            obj.plotting = plotting;
            if plotting
                drawplotter(obj), 
            else
                drawwaitbar(obj,0,'Initializing...','name',obj.name);
            end

            function C = getcolors()
            % Return persistent random color variations for sky/gnd/solar regions

                C.sky = randomcolormap([0.5 0.7 1],16,0.2);
                C.gnd = randomcolormap([0.6 0.8 0],8,0.2);
                C.sun = randomcolormap([1 1 0],4,0.2);
                C.sun(:,4) = 0.5;
                C.trck = zeros(0,3); % updated at runtime

                function c = randomcolormap(center,n,var)
                % Return n random color variations around center color
                    c = rand(n,3)*var - var/2;
                    c = center + c - mean(c,1);
                    c = min(max(c,0),1);
                end
            end
        end

        function drawplotter(obj,varargin)

            args = {'Shading-Analysis Plotter'};
            if obj.diffuseshading, args{2} = 'big'; else, args{2} = '2:1'; end
            drawplotter@plotter(obj,args{:});

            obj.subplothandles = struct();
            if obj.diffuseshading
                tags = {'scene','beamsh','shregions','povproj'};
                for j = 1:4
                    obj.subplothandles.(tags{j}) = subplot(2,2,j);
                    set(gca,'visible','off');
                    set(gca,'position',[0.06+(mod(j,2)==0)*0.46,0.06+(j<=2)*0.46,0.42,0.42]);
                end 
            else
                obj.subplothandles.scene = subplot(1,2,1);
                set(gca,'position',[0.06,0.1,0.4,0.8]);
                obj.subplothandles.beamsh = subplot(1,2,2);
                set(gca,'position',[0.54,0.1,0.4,0.8]);
            end
            
            [~,obj.updatefcn] = plottrackerarray(obj.Mounts,obj.sunaz,obj.sunel,obj.Gnd,...
                obj.groundshading,'mask',true,'btn',0,obj.subplothandles.scene);
        end

		function setaxislimits(obj,tr)
        % Determine axis limits from a set of tracker centers
            assert(ishandle(obj.fighandle),'shadeanalysisplotter:noobj','Closed by user');
            
            [lims(:,1),lims(:,2)] = bounds([obj.Mounts.centers(:,obj.Mounts.masks(tr,:)),...
                                            obj.Mounts.centers(:,tr)],2);
			scale = max(diff(lims,1,2)); % get largest system dimension
            % scale = max(scale, sqrt(area(obj.Mounts.geom.border))*10);
			origin = obj.Mounts.centers(:,obj.Mounts.analysedtrackers(tr,:));
            
			lims = repmat(origin,1,2)+[-2,2;-2,2;-1,1]*scale*0.3;
			lims = lims';
			obj.axislimits = lims(:);
		end

		function minorupdate(obj,t,tr,offset)
            assert(ishandle(obj.fighandle),'shadeanalysisplotter:noobj','Closed by user');
            validateattributes(offset,'numeric',{'finite','real','size',[3,1]});
            if ~isequal(obj.axisoffset,offset)
                obj.updatefcn(t,tr,offset);
            end
        end

		function majorupdate(obj,t,tr,BLPoly,STRI,BTRI)
        % Update layout geometry (and ground shades) at each timestep,...
        % Move axes and plot beam shades for each new mount
        
            if nargin < 5, BLPoly = []; end
            if nargin < 6, STRI = []; end
            if nargin < 7, BTRI = []; end
                    
            obj.setaxislimits(tr);
            obj.updatefcn(t,tr,[],isfield(obj.Mounts,'masks'),STRI,BTRI);
            % axis(obj.subplothandles.scene,obj.axislimits);
            
            title(obj.subplothandles.scene,...
                sprintf('sunel %0.1f^o, sunaz %0.1f^o (%s)',obj.sunel(t),obj.sunaz(t)));
            
            % Plot beam shades (if any)
            obj.plotbeamshades(BLPoly);
            title(obj.subplothandles.beamsh, ...
                sprintf('t = %03d, tr = %02d',t,tr));
        end

% function plotscenario(obj,TckP,mask,GndShP)
%             persistent ax;
%             persistent gshh;
%             persistent polyh;
%             
%             assert(ishandle(obj.fighandle),'shadeanalysisplotter:noobj','Closed by user');
%             figure(obj.fighandle);
% 			subplot(obj.subplothandles.scene); % cla();
%             
%             m = nnz(mask);
%             c = jet(m);
%             obj.colors.trck = c; % remember colors for plotprojectedtrck
%             colormap([c;[0.5,0.5,0.5]]);
%             c = zeros(numel(TckP),1);
%             c(~mask) = m+1; % gray
%             c(mask) = 1:m;
% 
% 			axisscale = max(diff(reshape(obj.axislimits,2,3),1,1))*0.15; 
% 
% 			hold on
% 			% Plot absolute axis (polyline: x 0 y 0 z)
%             if ~isempty(ax) && ishandle(ax), delete(ax); end
% 			ax = [1 0 0 0 0;0 0 1 0 0;0 0 0 0 1]*axisscale+repmat(obj.axisoffset,1,5);
% 			ax = plot3(ax(1,:),ax(2,:),ax(3,:),'k-');
% 
% 			% Plot Shade Polygons
%             if ~isempty(gshh) && any(ishandle(gshh)), delete(gshh); end
%             gshh = polyplot(GndShP,[0.8,0.8,0.8],'none');
% 
% 			% Plot Tracker Polygons
%             if isempty(polyh) || ~ishandle(polyh)
%                 [V,~,F] = poly2vef(TckP,true);
%                 polyh = patch('Faces',F,'Vertices',V,'FaceVertexCData',c,...
%                     'facecolor','flat','CDataMapping','direct');
%             else
%                 set(polyh,'Vertices',poly2vef(TckP));
%                 set(polyh,'FaceVertexCData',c);
%             end
% 			
% 			xlabel('x'); ylabel('y')
% 			axis equal; 
% 			view(120,20)
% 			axis(obj.axislimits)
% 			
% 				%axis(axis);
% 				% % Plot rough ground mesh:
% 				% [Tri,V] = groundmesh(Trackers);
% 				% colormap('summer');
% 				% trisurf(Tri,V(:,1),V(:,2),V(:,3));
% 			hold off
% 		end

		function plotbeamshades(obj,BLPoly)
        % Plot frontal view of (shaded) tracker under analysis
        
            MOUNT = [0.1 0.1 0.4];
            MOUNT_EDGE = 'none';
            LIGHT = [0.2 0.2 0.8];
            LIGHT_EDGE = [0.4 0.4 0.9];
            SHADE = [1 1 0 0.2];
            SHADE_EDGE = 'y';
        
            assert(ishandle(obj.fighandle),'shadeanalysisplotter:noobj','Closed by user');
            figure(obj.fighandle);
			subplot(obj.subplothandles.beamsh); 
			cla();
            trackeroutline = obj.Mounts.geom.border;
            if isstruct(trackeroutline), trackeroutline = polygon.unpack(trackeroutline); end

            if nargin < 2 || isempty(BLPoly) || isequal(BLPoly,0)
                polyplot(trackeroutline,MOUNT,MOUNT_EDGE);
            elseif isequal(BLPoly,1)
                polyplot(trackeroutline,LIGHT,LIGHT_EDGE);
                polyplot(trackeroutline,SHADE,SHADE_EDGE);
            elseif isstruct(BLPoly)
                polyplot(trackeroutline,MOUNT,MOUNT_EDGE);
                pvArea.shadingfactors(obj.Mounts.geom,BLPoly,'-pack16','depth',3,...
                    '-plot','alpha',1,'color',LIGHT,'edgecolor',LIGHT_EDGE);
                polyplot(polygon.unpack(BLPoly),SHADE,SHADE_EDGE);
                % p = polygon.unpack(PolyShadesPOA);
				% polyplot(p,[0.2 0.2 0.8],'none');
            else
                error('Unexpected PolyShadesPOA');
            end

			xlim = max(abs(trackeroutline.x));
			ylim = max(abs(trackeroutline.y));
			xlim = ceil(xlim + 0.2*min(xlim,ylim)); % add margin
			ylim = ceil(ylim + 0.2*min(xlim,ylim));

			axis equal
			axis([-xlim,xlim,-ylim,ylim])
		end

		function plotskyregions(obj,prjreg,shaded,W_sky,W_gnd,W_cs)
		% Plot proyections of any POV-independent elements of the scenario, 
        % i.e. horizon lines & sky regions
            
            for f = {'sky','albedo','solar';'sky','gnd','sun'}
                if numel(prjreg.(f{1})) > size(obj.colors.(f{2}))
                    cidx = mod(0:numel(prjreg.(f{1}))-1,size(obj.colors.(f{2}),1))+1;
                    obj.colors.(f{2}) = obj.colors.(f{2})(cidx,:);
                end
            end
            
            assert(ishandle(obj.fighandle),'shadeanalysisplotter:noobj','Closed by user');
            figure(obj.fighandle);
            
            if ~shaded
                subplot(obj.subplothandles.shregions);
            else
                subplot(obj.subplothandles.povproj); 
            end
            
            % Unpack polygons
            prjreg.sky = cellfun(@polygon.unpack,prjreg.sky,'unif',0);
            prjreg.albedo = cellfun(@polygon.unpack,prjreg.albedo,'unif',0);
            prjreg.solar = cellfun(@polygon.unpack,prjreg.solar,'unif',0);
            
            cla();
            polyplot(polygon(180),[1 1 1]*0.5,'none'); % black unit circle as background
            batchplot(prjreg.sky,obj.colors.sky);
            batchplot(prjreg.albedo,obj.colors.gnd);
            batchplot(prjreg.solar,obj.colors.sun);       
			% axis([-1.2 1.2 -1.2 1.2])
			axis square
            set(gca,'yticklabel',[])
			set(gca,'xticklabel',[])
            
            if nargin > 3
            % Text labels for DwF or DshF
                UINT16MAX = 2^16-1;
                lbl = @(F) arrayfun(@(n,f) sprintf('%d\n%0.2f',n,f),1:numel(F),double(F(:)')/UINT16MAX,'unif',0);
                
                x = labelpolygons(prjreg.solar,lbl(W_cs),[],mean(obj.colors.sun,1)*0.6);
                x = labelpolygons(prjreg.sky,lbl(W_sky),x,mean(obj.colors.sky,1)*0.6);
                    labelpolygons(prjreg.albedo,lbl(W_gnd),x,mean(obj.colors.gnd,1)*0.6);           
            end
            
            function batchplot(P,C)
            % Plot polygons P as patches with colors C
            
                if isempty(P), return; end

                haveholes = cellfun(@(p) any([p.hole]),P);
                if any(haveholes)
                    warning_reseter = naptime('MATLAB:delaunayTriangulation:ConsSplitPtWarnId'); %#ok<NASGU>
                    arrayfun(@(k) polyplot(P{k},C(k,:),'w'),find(haveholes));
                    batchplot(P(~haveholes),C(~haveholes,:));
                    return;
                end

                if size(C,2) == 4
                    alpha = mean(C(:,4));
                    C = C(:,1:3);
                else
                    alpha = 1;
                end

                idx = repelem(1:numel(P),cellfun(@numel,P));
                [V,~,F] = poly2vef([P{:}],'unif');
                patch('Faces',F,'Vertices',V,'FaceVertexCData',C(idx,:),'EdgeColor','w',...
                        'FaceColor','flat','FaceAlpha',alpha,'LineWidth',0.1);
            end
            
            function x = labelpolygons(P,L,x0,C)
            % Find maximum-inscribed-circle positions to place labels for a set of polygons

                if nargin < 2 || isempty(x0), x0 = zeros(0,2); end
                x = x0;
                for k = 1:numel(P)
                    if isvoid(P{k}), continue; end
                    c = maxinscribedcircle(P{k},x);
                    if ~isempty(c)
                        x = cat(1,x,c);
                        text(c(1),c(2),L{k},'color',C,'fontsize',8,...
                            'VerticalAlignment','middle','HorizontalAlignment','center');
                    end
                end
            end
        end

% 		function plotprojectedtrck(obj,PjTrckP,PjGnd  ShP)
% 		% Plot proyections of the trackers and ground shades, after removing any projections 
% 		% plotted on past iterations.
% 
%             assert(ishandle(obj.fighandle),'shadeanalysisplotter:noobj','Closed by user');
%             figure(obj.fighandle);
% 			subplot(obj.subplothandles.shregions);
% 			hold all
%             notempty = ~cellfun(@isempty,PjTrckP);
%             set(gca,'colororder',obj.colors.trck(notempty,:));
% 
% 			% Delete any plotted shades/trackers from past iterations
%             arrayfun(@deleteifhandle,obj.polyh.pjshades);
% 			obj.polyh.pjshades = NaN;
%             arrayfun(@deleteifhandle,obj.polyh.pjtrackers);
% 			obj.polyh.pjtrackers = NaN;
% 			
% 			% Plot the new set of shade/tracker projections
%             pplot = @(p,varargin) plot(p.x([1:end,1]),p.y([1:end,1]),varargin{:});
%             obj.polyh.pjtrackers = cellfun(pplot,PjTrckP(notempty));
%             obj.polyh.pjshades = arrayfun(@(x) pplot(x,'k-'),PjGndShP);
% 			hold off
%             
%             function deleteifhandle(x)
%                if ishandle(x), delete(x); end 
%             end
%         end
		
		function exportfigures(obj,t,tr,pt)
        % Save all plotted axes as .png figures, inside a ../Frames directory
            assert(ishandle(obj.fighandle),'shadeanalysisplotter:noobj','Closed by user');
            figure(obj.fighandle);
			ID = sprintf('%02d%d%03d.png',tr,pt,t);
			specs = {'-r180','-png','-transparent'};
			export_fig(sprintf('./Frames/scen%s',ID),specs{:},obj.subplothandles.scene);
            export_fig(sprintf('./Frames/beam%s',ID),specs{:},obj.subplothandles.beamsh);
            if obj.diffuseshading
                export_fig(sprintf('./Frames/proj%s',ID),specs{:},obj.subplothandles.povproj);
                export_fig(sprintf('./Frames/poly%s',ID),specs{:},obj.subplothandles.shregions);
            end
        end
    end
end
