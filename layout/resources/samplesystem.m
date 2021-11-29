function [Trck, arraydef] = samplesystem(type,varargin)
% [Trck, arrdf] = SAMPLESYSTEM(TYPE)
% [Trck, arrdf] = SAMPLESYSTEM(TYPE,MDIMS)
% [Trck, arrdf] = SAMPLESYSTEM(TYPE,MDIMS,PHYSPARAMS,ELECPARAMS)
% [Trck, arrdf] = SAMPLESYSTEM(..,NAME,VALUE,..)
%
% Generates example PV-Plant-layouts (and array-definitions) for use in simulation test runs.
%
% TYPE: string {'2a','1aV','1aF','1aC','0a'} (xa = n-axis, V = vertical, F = Fix-tilt, C = Contour)
% MDIMS: module dimensions [w,h] or polygon/pvArea object. Default calls stdmodulepoly().
% PHYSPARAMS: vector [n, m, N, M] - N rows of M mounts (each n rows of m modules) 
% ELECPARAMS: vector [sn,sm,sN,sM] - the [n·N,m·M] array of modules is connected on strings of 
%       sn rows and sm modules. Each sN x sM array of strings is connected to a different inverter.
%       Verify that n*N/sN and m*M/(sM*sm) are integers. Use {} to skip electrical definition.
%
%     Current defaults yield a 1600 module array, roughly 80 m x 80 m in size, connected in 
%     20-module strings, 20-strings-per-inverter, for all mount types.
%
% SAMPLESYSTEM(...,'pitch',p) use pitch p between rows. Default for p is twice table height.
% SAMPLESYSTEM(...,'pitch',[p,q]) use pitch p between rows, q between columns.
%       Default for q is table width plus tgap between mounts (see below)
%
% SAMPLESYSTEM(...,'az',az) rotate array by azimuth az. Default is zero.
%
% SAMPLESYSTEM(...,'tilt',th) use fixed tilt th for {0a,1aV,1aF} tracker types. Defaults between
%   0-30° will be set depending on mount type.
%
% SAMPLESYSTEM(...,'mgap',gap) use scalar or 2-vector gap [gx,gy] between modules.
%       Default is 2% of shortest module edge.
%
% SAMPLESYSTEM(...,'tgap',gap) use (scalar) gap between mounts when q is not provided as second
%       element of 'pitch' option. Default is 20% of shortest module edge.
%
% SAMPLESYSTEM(...,'location',Loc) use Loc.longitude, Loc.latitude as layout-origin
%       otherwise the system will be placed on Null Island (0°, 0°)
%
% SAMPLESYSTEM(...,'landscape',f) where f is boolean (false = flat, true = sample 3D) or a 
%       function handle f(x,y). Use a surface function to generate 3D landscape.
%
% SAMPLESYSTEM(...,'analysedpts',[c,r]) set TRCK.analysedpoints = QUADRATUREPOINTS([c,r],A,B)
%        where A, B are the mount's rectangular limits. Default is [-0.4,-1.1]
%                                        
% SAMPLESYSTEM(...,'files',name) generates pre-processor-output-like files for the sample system,
%       According to Pre-Processor Specifications V2.0. File names follow the conventions:
%
%       name_980x1960mm_72c_3bd.mpoly - Module file: 72-cell, 3 bypass-diodes, 0.98 x 1.96 m
%       name_1R20VM_980x1960mm.tpoly - Mount file: 1-row of 20-Vertical-modules per mount
%       name_20R4T_1R20VM_980x1960mm.mounts - Layout file: 20-rows of 4 mounts
%       name_20R4T_1R20VM_8i10s20M.arrdef - Array Def.: 8-Inverters, each 10 strings of 20 modules
%
%   Additional configuration details (pitch, azimuth, gaps, etc.) included in header-comments.
%
% SAMPLESYSTEM(...,'files',true) uses default-generated name.
%
% See also: STDMODULEPOLY, SAMPLEARRAYDEF, WRITETRACKERLAYOUT, EXPORTPOLYGONSFILE,
%           WRITEARRAYDEFINITION, QUADRATUREPOINTS, TESTBASE

%% Parse input and options

    opt.mdims = [];
    opt.physparams = [];
    opt.elecparams = [];
    opt.location = struct('latitude',0,'longitude',0,'altitude',0,'TimeZone',0); 
    opt.rotation = [];
    opt.pitch = [];
    opt.az = 0;
    opt.tilt = [];
    opt.mgap = [];
    opt.tgap = [];
    opt.files = false;
    opt.landscape = false;
    opt.analysedpts = [-0.4000 -1.1000];
    opt.gndclearance = 0.5;
    
    % Check the remaining arguments for options, use defaults if missing (*)
    [opt,~,isdef] = getpairedoptions(varargin,opt,'dealrest',4);

    [mdims,physparams,elecparams] = deal(opt.mdims,opt.physparams,opt.elecparams); 
    if ~isdef.location, opt.location = parselocation(opt.location,'-soft'); end

    % Check type, and fix any upper/lower-case issues
    typelist = {'2a','1aV','1aF','1aC','0a'};
    type = parselist(type,typelist,'tracker type');
    
    Trck.name = 'foo'; % just to make it the first field
    Trck.type = type;
    Trck.origin = [opt.location.longitude,opt.location.latitude,opt.location.altitude];
    if isempty(opt.rotation)
        Trck.rotation = opt.az + 180*(opt.location.latitude < 0);
        opt.az = 0;
    else
        validateattributes(opt.rotation,{'numeric'},{'real','scalar','>=',-360,'<=',360});
        Trck.rotation = opt.rotation;
    end
    
    if isempty(opt.tilt)
        switch type
            case '1aV', opt.tilt = max(20,abs(opt.location.latitude)); % 1aV
            case {'1aC','1aF'}, opt.tilt = min(30,abs(opt.location.latitude)); % 1aC / 1aF
            case '0a', opt.tilt = max(10,min(30,abs(opt.location.latitude))); % 0a
        end
    end

    % Get module dimensions and 3-depth pvArea from whatever is in MDIMS
    if isempty(mdims)
        modarea = stdmodulepoly();
        [mw,mh,~] = rectangleproperties(modarea.border); mdims = [mw;mh];
        info.module = 'STDMODULEPOLY()';
    elseif isa(mdims,'polygon')
        [mw,mh,~] = rectangleproperties(mdims); mdims = [mw;mh];
        modarea = [];
        info.module = 'POLYGON object';
    elseif isa(mdims,'pvArea')
        modarea = mdims;
        [mw,mh,~] = rectangleproperties(modarea.border); mdims = [mw;mh];
        info.module = 'PVAREA object';
    elseif isnumeric(mdims) && numel(mdims) == 2
        mdims = mdims(:);
        modarea = [];
        info.module = 'MDIMS';
    else
        error('sampletrackers:mdims','Expecting 2-vector or pvArea as second argument');
    end
    if isempty(modarea)
        try modarea = stdmodulepoly([],mdims); 
        catch ERR
            error(['STDMODULEPOLY failed to find a standard cell-layout for a module with ',...
                   'dimensions %0.1f×%0.1f mm, generate a PVAREA object first and pass it as ',...
                   'an argument to SAMPLESYSTEM. ERR: %s'],mdims*1000,getReport(ERR));
        end
    end

    % Replace defaults by actual parameters, if available
    if isempty(physparams)
        % Default array dimensions/characteristics. 
        % Current numbers yield a 1600 module array, roughly 80 m x 80 m in size, connected in 
        % 20-module strings, 20-strings-per-inverter, for all mount types.
        % If you plan on changing anything, check that n·N/(sN·sn) and m·M/(sM·sm) are integers, 
        % and consider using constant n·m·N·M, sm·sn, and sN·sM to use the same inverter
        switch type
            case '0a'
                physparams = [4  5 10 8]; % [n, m, N, M] - H, 10x4 m tracker
                elecparams = [1 20 20 1]; % [sn sm sN sM] 
            case {'1aV','2a'}
                physparams = [8 5 5 8]; % [n, m, N, M] - H, 10x8 m tracker
                elecparams = [4 5 5 4]; % [sn sm sN sM] 
            case '1aF'
                physparams = [2 10 10 8]; % [n, m, N, M] - V, 10x4 m tracker
                elecparams = [1 20 10 2]; % [sn sm sN sM] 
            case '1aC'
                physparams = [1 20 20 4]; % [n, m, N, M] - V, 20x2 m tracker
                elecparams = [1 20 10 2]; % [sn sm sN sM] 
        end
        
        % When using default configurations, adjust module orientation as well
        switch type
            case {'1aC','1aF'}, ok = mdims(1) <= mdims(2);    % Vertical module (h > w)
            otherwise, ok = mdims(1) >= mdims(2);             % Horizontal (w > h)
        end
        if ~ok
            modarea = rotate(modarea,90);
            mdims = mdims([2,1]);
        end
    end
    args = num2cell(physparams);
    [n, m, N, M] = deal(args{:});
    
    % Get module short-label and description
    [modulelbl,info.module] = moduleinfo(mdims,modarea,info.module);
    
    if isempty(opt.mgap), opt.mgap = round(min(mdims)*0.02,3); end   %(*) complete defaults
    if isempty(opt.tgap), opt.tgap = round(min(mdims)*0.20,2); end                           
    if isscalar(opt.mgap), opt.mgap = [1;1]*opt.mgap; else, opt.mgap = opt.mgap(:); end
    
    [mountlbl,info.mount] = mountinfo(n,m,mdims,type);
    [layoutlbl,info.layout] = layoutinfo(N,M,mountlbl,info.mount);
        
%% Build the actual system
    
    % Get table dimensions
    trckdims = [m;n].*mdims + ([m;n]-1).*opt.mgap;   % mount size

    % Use default pitch and horizontal-spacing, if necessary
    % Defaults can be type-specific, considering typelist = {'2a','1aV','1aF','1aC','0a'}
    if ~isempty(opt.pitch), p = opt.pitch(1);
    else
        p = round(strcmp(type,typelist)*[3 3 2.5 2.5 2]'*trckdims(2),1); 
        opt.pitch = p;
    end
    if ~isscalar(opt.pitch), q = opt.pitch(2);
    else
        q = round(strcmp(type,typelist)*([2 2 1.5 1 1]*trckdims(1) + [0 0 0 1 1]*opt.tgap)',1);
        opt.pitch(2) = q;
    end

    % Generate array definition, if required
    if ~isempty(elecparams)
        elecparams = num2cell(elecparams);
        [sn,sm,sN,sM] = deal(elecparams{:});
        arraydef = samplearraydef(n,m,N,M,sn,sm,sN,sM);  
        arraydefexists = true;
        [arraydeflbl,info.arraydef] = arraydefinfo(n,m,N,M,sn,sm,sN,sM);
    else
        arraydefexists = false;
    end
        
    % Create pvArea representing table/tracker
    % a = pvArea(mdims(1),mdims(2));  % simplified unit element                 
    a = modarea;
    step = opt.mgap + mdims; 
    cx = mdims(1)/2+(0:m-1)*step(1)-trckdims(1)/2;  % module-center vectors
    cy = mdims(2)/2+(0:n-1)*step(2)-trckdims(2)/2;
    
    sa = areaArray(a,cx,zeros(size(cx)));
    Trck.geom = areaArray(sa,zeros(size(cy)),cy);
    %[cy,cx] = meshgrid(cy,cx);        % module numbers along module rows (left->right, down->top)
    
    %DEBUG: add "arrows" at x', y'
    % Trck.geom.border = makepointy(Trck.geom.border);

    % Generate array of tracker centers
    [Y,X] = meshgrid(-p/2:p:(N-1)*p,-q/2:q:(M-1)*q);  % arrange mount numbers along rows
    X = X(:); Y = Y(:);
    X = X - mean(X); Y = Y - mean(Y);   % Center around origin

    if any(strcmp(type,{'1aC','1aF'}))
        [X,Y] = deal(-Y,X); % rotate -90°
    end
        
    Trck.centers = [X,Y,zeros(size(X))]'; % Trck.centers stay in (rotated) project coordinates

    % Rotate by r to get real coordinates
    r = Trck.rotation;
    [X,Y] = deal(X*cosd(r)-Y*sind(r),X*sind(r)+Y*cosd(r));

    if isequal(opt.landscape,true)
    % generate a sample landscape, if requested
        opt.landscape = @(x,y) 4*cos(sqrt((3*x/60-2).^2 + (2*y/60+1).^2));
    end
    
    % GUIfigure('samplesystem'); clf(); hold on;
    % [xx,yy] = meshgrid(-50:50,-50:50);
    % contour(xx,yy,opt.landscape(xx,yy))
    % plot(X,Y,'r.')
    % arrayfun(@(x,y,n) text(x,y,num2str(n)),X,Y,(1:numel(X))')
    % axis equal;

    if isequal(opt.landscape,false)
    % everything is flat
        axistilt = zeros(size(X));
    else
    % try calling landscape as a function, to accept interpolants
        try assert(isnumeric(opt.landscape(0,0)));
        catch, error('bad landscape: expecting logical value or function-handle');
        end
        switch type
            case {'0a','1aC'}
                % centerpoint is the average of extemes, not the function evaluated at center
                r = Trck.rotation; % x' axis
                if isequal(type,'1aC'), r = r + 90; end
                q = q*[cosd(r),sind(r)];
                Z(:,1) = opt.landscape(X+q(1)/2,Y+q(2)/2);
                Z(:,2) = opt.landscape(X-q(1)/2,Y-q(2)/2);
                axistilt = atand((Z(:,1)-Z(:,2))/norm(q));
                Zc = (Z(:,1)+Z(:,2))/2;
            otherwise
                Zc = opt.landscape(X,Y);
        end
        Trck.centers(3,:) = Zc';
    end

    % Common propperties for 1aF and 1aC
    if any(strcmp(type,{'1aF','1aC'}))
        Trck.azimuth = opt.az;
        Trck.axisoffset = [0,0,0.1];
        Trck.tracklimits = [-45,45];
        Trck.groundcoverratio = trckdims(2)/p;
        Trck.backtracking = true;
        Trck.group = reshape(ones(M,1)*(1:N),[],1); % group by rows
    end

    % Assign type-specific propperties
    switch type
        case '2a'
            Trck.tracklimits = [-150,150,5,90];
            Trck.axisoffset = [0,0,0.5];
            opt = rmfield(opt,'tilt'); % not used
        case '1aV'
            Trck.tilt = opt.tilt;
            Trck.tracklimits = [-150,150];
            Trck.axisoffset = [0,0,0.5];
        case '1aF'
            Trck.tilt = opt.tilt;
        case '1aC'
            Trck.tilt = axistilt;
            Trck.groupidx = reshape(repmat(1:M,N,1),1,[]); % set rows as groups
            opt = rmfield(opt,'tilt'); % not used
        case '0a'
            Trck.axisoffset = [0,0,0];
            Trck.slope = axistilt;
            Trck.tilt = opt.tilt;
            Trck.azimuth = opt.az;
    end
    
    Trck.centerheight = eps();
    Trck.centerheight = opt.gndclearance - min(checkmountclearance(Trck,-Inf));
    Trck.centerheight = ceil(Trck.centerheight*20)/20; % round to 5 cm
    
    % Define diffuse-shading-analysis nodes
    Trck.analysedpoints = round(quadraturepoints(opt.analysedpts,-trckdims/2,trckdims/2)',2);
        
    % % analyze up to NaxNa representative trackers, distributed as Legendre-Gauss points
    % Na = 4;
    % repidx = quadraturepoints(Na,[1,1],[M,N]);
    % repidx = unique(round(repidx),'rows');
    % analysed = false(M,N); 
    % analysed(repidx(:,1),repidx(:,2)) = true;
    % Trck.analysedtrackers = find(analysed);
 
    % Generate names and labels
    if ischar(opt.files)
        prefix = opt.files;
        opt.files = true;
    else, prefix = 'sample'; 
    end
    
    Trck.info = optioncellstr(opt);
    Trck.info = cat(1,{info.layout},Trck.info(:));
            
    mountlbl = [prefix '_' mountlbl];
    layoutlbl = [prefix '_' layoutlbl];
    if arraydefexists, arraydeflbl = [prefix '_' arraydeflbl]; end

    Trck.name = layoutlbl;
    Trck.modelidx = LabelMap(1,{mountlbl});     % PROVISIONAL: single mount model

%% Export files, if requested

    if opt.files
        
        listvector = @(c) regexprep(mat2str(c(:)),'[[]]','');
                
        modspecs.module_id = modulelbl;
        modspecs.comments = {info.module};
        filenames{1} = [modulelbl '.mpoly'];
        
        trckspecs.mount_model_id = mountlbl;
        trckspecs.module_id = modulelbl;
        trckspecs.mount_type = type;
        trckspecs.axis_offset = listvector(Trck.axisoffset);
        trckspecs.center_height = num2str(Trck.centerheight);
        trckspecs.analysed_pts = listvector(Trck.analysedpoints);
        if isfield(Trck,'tracklimits'), trckspecs.track_limits = listvector(Trck.tracklimits); end
        if isfield(Trck,'backtracking'), trckspecs.backtracking = num2str(Trck.backtracking); end
        if isfield(Trck,'tilt') && isscalar(Trck.tilt), trckspecs.tilt = num2str(Trck.tilt); end
        
        filenames{2} = [mountlbl '.tpoly'];
        filenames{3} = [layoutlbl '.mounts'];

        if arraydefexists, filenames{4} = [arraydeflbl '.arrdef']; end
        
        % Check if files exist, and if so, if they can be overwritten
        ok2write = right2overwrite(filenames); 
        if ~ok2write; return; end

        % Write module file
        exportpolygonsfile(modarea,filenames{1}, modspecs, ok2write);
        
        % Write layout file
        writetrackerlayout(Trck,ok2write);

        % Write mounts file
        trckspecs.comments = {'Generated with SAMPLESYSTEM(...)'; info.mount; Trck.info{2}};
        exportpolygonsfile(reducedetail(Trck.geom,3),filenames{2}, trckspecs, ok2write);
        
        % Write array definition
        if arraydefexists
            arrdefspecs = struct('arrdef_id',arraydeflbl,'layout_id',layoutlbl);
            arrdefspecs.comments = {'Generated with SAMPLESYSTEM(...)'; info.arraydef};
            writearraydefinition(arraydef,filenames{4},arrdefspecs,ok2write); 
        end
    end
end

function s = optioncellstr(opt)
% List layout options, to include as comment in export-files

    s{1} = sprintf('module-gaps = %0.1f × %0.1f cm',opt.mgap*100);
    s{2} = sprintf('mount-gap = %0.1f cm, pitch = %0.1f × %0.1f m',opt.tgap*100,opt.pitch);
    s{3} = sprintf('azimuth = %0.1f°', opt.az);
    if isfield(opt,'tilt'), s{3} = sprintf('%s, tilt = %0.1f°', s{3}, opt.tilt); end
    if isequal(opt.landscape,0), s{4} = 'flat landscape';
    elseif isa(opt.landscape,'function_handle')
        s{4} = ['landscape function = ' func2str(opt.landscape)];
    else
        s{4} = ['landscape from' class(opt.landscape)];
    end
end

function [modulelbl,info] = moduleinfo(mdims,modarea,source)
% Generate module short-label and description
% modulelbl: '980x1960mm_72c_3bd'
% info: '72-cell (3-bypass-diode) 980x1960mm module from SOURCE'

    nc = prod(modarea.dims);
    dimlbl = sprintf('%dx%dmm',round(mdims(1)*1000),round(mdims(2)*1000)); % e.g. 980x1960mm
    if modarea.depth < 3, nd = 0; else, nd = modarea.dims(1); end
    modulelbl = sprintf('%s_%dc_%dbd',dimlbl,nc,nd);
    info = sprintf('%d-cell (%d-bypass-diode) %s module from %s',nc,nd,dimlbl,source);
end
    
function [mountlbl,info] = mountinfo(n,m,mdims,type)
% Generate mount short-label and description
% mountlbl: '0a_1R20VM_980x1960mm'
% info: 'Fixed tables, 1 row of 20 vertical modules (980x1960mm)'

    if mdims(1) < mdims(2), side = 'vertical'; else, side = 'horizontal'; end
    
    mountlbl = sprintf('%dR%d%cM',n,m,upper(side(1)));                       % 1R20VM
    dimlbl = sprintf('%dx%dmm',round(mdims(1)*1000),round(mdims(2)*1000));   % 980x1960mm
    mountlbl = strjoin({type,mountlbl,dimlbl},'_');                          % 0a_1R20VM_980x1960mm
    
    switch type
        case '0a' , typestr = 'Fixed tables';
        case '1aC', typestr = 'Horizontal(contour)-axis trackers';
        case '1aF', typestr = 'Fixed-tilt single-axis trackers';
        case '1aV', typestr = 'Vertical-axis trackers';
        case '2a' , typestr = 'Two-axis trackers';
    end
    info = sprintf('%s, %s of %s (%s)',typestr,nthings(n,'row'),nthings(m,[side ' module']),dimlbl);
end

function [layoutlbl,info] = layoutinfo(N,M,mountlbl,mountinfo)
% Generate layout short-label and description
% layoutlbl: '20R4T_0a_1R20VM_980x1960mm'
% info: '20 rows of 4 mounts: MOUNTINFO'

    layoutlbl = sprintf('%dR%dT_%s',N,M,mountlbl);
    info = sprintf('%s of %s: %s',nthings(N,'row'),nthings(M,'mount'),mountinfo);
end
    
function [arraydeflbl,info] = arraydefinfo(n,m,N,M,sn,sm,sN,sM)
% Generate array-definition short-label and description
% arraydeflbl: 2x2i_20x1s_1x20M
% info: '2x2 inverter-blocks with 20 strings of 20 modules (1 row)'

    arraydeflbl = sprintf('%dx%di_%dx%ds_%dx%dM',N*n/(sN*sn),M*m/(sM*sm),N*n/sn,M*m/sm,sn,sm);
    
    % Choose inverter-type
    if sn*sm*sN*sM == 1, info = 'micro-inverter';
    elseif sN*sM == 1, info = 'string-inverter'; else, info = 'inverter-block'; 
    end
    
    % Fit inverter number
    if n*m*N*M == sn*sm*sN*sM, info = ['1 ' info];
    else, info = sprintf('%dx%d %ss',n*N/(sN*sn), m*M/(sM*sm),info);
    end
    
    % Add N° of strings, for multi-string inverters
    if sN*sM > 1, info = [info ' with ' nthings(sN*sM,'string')]; end

    % Add N° of modules, unless micro-inverters are used
    if sn*sm*sN*sM > 1
        info = sprintf('%s of %d modules (%s)',info,sn*sm,nthings(sn,'row'));
    end
end

% function P = makepointy(P)
%     [w,h] = rectangleproperties(P);
%     t = polygon(3);
%     xt = polytranslate(t,[w/2,0]);
%     yt = polytranslate(polyrotate(t,90),[0,h/2]);
%     P = mergepolygons([P,xt,yt],'pos');
% end