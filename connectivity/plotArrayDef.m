function varargout = plotArrayDef(ArrayDef,Mounts,varargin)
% PLOTARRAYDEF(ARRAYDEF,MOUNTS) - create a color-coded representation of the connection scheme
%   ARRAYDEF for the plant in MOUNTS.
% PLOTARRAYDEF(ARRAYDEF,MOUNTS,ARRIDX) - Use structure of NDmap objects ARRIDX for labeling of 
%   array elements. ARRIDX is expected to have fields ('i' and 's') or 'is', for inverter, string,
%   or inverter-string labels. 
% and optionally use figure H instead of the default figure-handle.

    opt.fh = [];
    opt.equator = 'up';
    opt.maxdetail = 3;
    opt.fontsize = [10,7];
    opt.labeling = []; % based on opt.maxlabels
    opt.maxlabels = 1000;
    opt.flip = false;
    [opt,remargs] = getpairedoptions(varargin,opt);
    assert(numel(remargs)<=1,'Unrecognized arguments');
    
    if ~isempty(remargs), ArrIdx = remargs{1}; else, ArrIdx = []; end
    
    centers = Mounts.centers(1:2,:);
    N = size(centers,2);
    
    R = mountrotations(Mounts,90,90);
    if size(R,3) == 1, R = repmat(R,[1,1,N]); end
	
    a0 = reducedetail(Mounts.geom,2);
    trckcopy = @(j) translate(rotate(copy(a0),flatrotation(R(:,:,j))),centers(:,j));
    TrckField = pvArea();			   
    TrckField.elements = arrayfun(trckcopy,1:N);
	TrckField = flattenpvarea(TrckField);
	
	connected = ismember(1:numel(TrckField.elements),ArrayDef.pidx(:));
    % pidx generates indices compatible with flattenpvarea - [Np,Nt]

    esubs = fliplr(ArrayDef.esubs(:,:));
    
    % % Remove trailing singleton dimensions (highest-hierarchy)
    % while size(esubs,2)>1 && all(esubs(:,end)==esubs(1,end))
    %     esubs = esubs(:,1:end-1);
    %     opt.maxdetail = opt.maxdetail-1;
    % end
	ConnTrck = reshapearea(TrckField.elements(connected),esubs);
    ConnTrck.border = polygon();
	% writearraydefinition(ConnTrck,'arraydef.xlsx');
    
    if isempty(opt.fh), opt.fh = GUIfigure('arraydef','Array Definition','big'); end
    if nargout > 0, varargout{1} = opt.fh; end
    
	figure(opt.fh); clf;
    set(gcf,'Name','Array Definition'); set(gcf,'Numbertitle','off');
    switch lower(opt.equator)
        case {'up','+y'}
        case {'down','-y'}
            set(gca,'Xdir','reverse') % Equator (+y) facing down
            set(gca,'Ydir','reverse')
        otherwise, error('Unknown ''equator'' option value');
    end
    
    plotelements(ConnTrck,opt.maxdetail+1);
    if isempty(opt.labeling)
        j = 0;
        while j < opt.maxdetail && prod(ConnTrck.dims(1:j+1)) < opt.maxlabels
            j = j+1;
        end
        labeling = j;
    else
        labeling = opt.labeling;
        assert(isnumeric(opt.labeling) && opt.labeling <= opt.maxdetail && mod(opt.labeling,1) == 0,...
            '''labeling'' must be an integer <= ''maxdetail''')
    end
	plotArea(ConnTrck,opt.maxdetail,isempty(ArrIdx) && labeling > 0,0.2);
    
    if opt.flip, set(gca,'xdir','reverse','ydir','reverse'); end
    set(gca,'visible','off');
    equalunits(); % axis equal, but not square

    if ~isempty(ArrIdx) && labeling > 0
        if ~isfield(ArrIdx,'is'),ArrIdx.is = LabelMap(compose(ArrIdx.i,ArrIdx.s)); end
        
        Ni = ArrayDef.esize(1);
        Ns = ArrayDef.esize(2);
        assert(all([Ni,Ns]==ArrIdx.is.esize),'plotArrayDef:arridxesize',...
            'ArrIdx.esize ~= ArrayDef.esize(1:2)');
        
        if labeling > 1
            P{1} = [ConnTrck.elements.border];
            lbl{1} = ArrIdx.is.plabels(:,1); 
            lbl{1} = lbl{1}(:,1);
            for k = Ni:-1:1
                P{k+1} = [ConnTrck.elements(k).elements.border];
                lbl{k+1} = ArrIdx.is.plabels(k,:);
                lbl{k+1} = lbl{k+1}(:,end);
            end
            P = cat(2,P{:})';
            lbl = cat(1,lbl{:});
        else
            P = [ConnTrck.elements.border];
            lbl = ArrIdx.is.plabels(:,1); lbl = lbl(:,1);
        end
        
        format = {'horizontalalignment','center','verticalalignment','middle',...
                  'fontsize',opt.fontsize(2),'color','k'}; 
        h = labelpolygons(P,lbl,format{:});
        set(h(1:Ni),'fontsize',opt.fontsize(1),'fontweight','bold');
    end
    
    if all(connected), return; end
    for j = find(~connected(:))'
		p = patch(TrckField.elements(j).border.x,TrckField.elements(j).border.y,[0.8,0.8,0.8]);
        %set(p,'FaceColor',[0.8,0.8,0.8])
        %set(p,'FaceAlpha',1);
        set(p,'LineWidth',0.1);
        set(p,'EdgeColor',[0.6,0.6,0.6]);
        %set(p,'EdgeAlpha',1);
    end
        
    % % Print
    % view(180,90)
    % rez = 600; 
    % figpos = getpixelposition(opt.fh);
    % resolution = get(0,'ScreenPixelsPerInch');
    % set(opt.fh,'paperunits','inches','papersize',figpos(3:4)/resolution);...
    % set(opt.fh,'paperposition',[0 0 figpos(3:4)/resolution]);
    % print(opt.fh,'arrdef.png','-dpng',['-r',num2str(rez)],'-opengl')
end

function Q = flatrotation(R)
    x = R(1:2,1); x = x/hypot(x(1),x(2));
    y = R(1:2,2);
    y = y - (y'*x)*x; y = y/hypot(y(1),y(2));
    Q = [x,y];
end

function equalunits()
% axis equal, without forcing axis to be square

    axis('tight');
    L = axis();
    c = (L*[1 1 0 0; 0 0 1 1]')/2;
    dL = L*[-1 1 0 0; 0 0 -1 1]'*1.05;
    
    set(gca,'units','pixels');
    S = get(gca,'position')*[0 0 1 0; 0 0 0 1]';
    scale = max(dL./S);
    dL = S*scale;
    L = c+dL.*[-1;1]/2;
    axis(L(:)');
end