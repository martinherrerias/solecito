function testbase(varargin)
% TESTBASE(..) - Generate and [partially] run project files for each tracker type in SAMPLESYSTEM()
%
% '-cleanup': remove *.mounts/*.arrdef/*.tpoly/*.mpoly files in TYPE subdirectories
% '-restart': run splitscript('base.ssp',{'setup','models','meteo'}) on ./common, then replace
%       *.ssp project files on TYPE subdirectories. Forces 'cleanup' = true.
% '-solve': run projects from start to end, otherwise stop before shading
% '-flat': flat ground instead of artificial landscape
% 'files2copy', C: copy files C from ./common into every TYPE subdirectory
% 'types',C: restrict testbase to some mount types
% 'onerror', key: for key in {'debug','stop','continue'} choose what to do if an error is caught
%   during execution.

    % Parse options
    [opt,varargin] = getflagoptions(varargin,{'-cleanup','-restart','-solve','-flat'});
    opt.files2copy = {};
    opt.types = {'0a','1aC','2a','1aV','1aF'};
    if ~isempty(dbstatus), opt.onerror = 'debug'; else, opt.onerror = 'stop'; end
    [opt,varargin] = getpairedoptions(varargin,opt,'restchk');
    
    if opt.restart, opt.cleanup = true; end
    if isempty(opt.files2copy), opt.files2copy = cell(0,1); end
    if ischar(opt.types), opt.types = {opt.types}; end
    
    switch lower(opt.onerror)
        case 'stop'
            stoponerror = true;
        case 'debug'
            dbstop if error;
            stoponerror = false;
        case 'continue'
            stoponerror = false;
        case 'exit' 
            stoponerror = false; % never mind, 'exit' goes to splitscript
        otherwise
            error('Unrecognized ''onerror'' option.');
    end
    
    if ~isempty(varargin)
        assert(numel(varargin) == 1 && ischar(varargin{1}),'');
        basepath = varargin{1};
    else
        basepath = pwd(); % fileparts(mfilename('fullpath'));
    end
    here = pwd(); lastwill = onCleanup(@() cd(here));
    cd(basepath);

    if opt.cleanup && any(cellfun(@isfolder,opt.types))
        if ~strcmp(optquestdlg(['TESTBASE will remove any *.mounts/*.arrdef/*.tpoly/*.mpoly files in ',...
                    shortliststr(opt.types(cellfun(@isfolder,opt.types)),'folder')],'Warning',...
                    'Go for it!','Wait, what? No!','Go for it!'),'Go for it!')
              return;
        end
    end
    
    fprintf('Running TESTBASE on %s\n\n',basepath);
    disp(opt);
    
    if ~isfolder('common'), mkdir('common'); end
    cd('./common');
    
    if opt.flat, optargs = {'flathorizon',true}; else, optargs = {}; end
    
    if opt.restart || isempty(dir('base.ssp'))
        if ~isempty(dir('base.log')), delete('base.log'); end
        f = splitscript('base.ssp',{'setup','models','meteo'},'onerror',opt.onerror,'simoptions',optargs);
        assert(~stoponerror || f == 0,'error while evaluatig testbase');
    end

    F = load('base.ssp','-mat','Location','ModIVint');
    mdims = [F.ModIVint.source.width, F.ModIVint.source.length]/1000;
    Loc = F.Location;
    
    if ~opt.flat
        %getgeotiff(Loc.latitude,Loc.longitude,50000);
        DEMfile = pickfile('*.tif');
        if isempty(dir(DEMfile)), copy(DEMfile,pwd); end
        DEM = getlandscape(DEMfile,Loc.latitude,Loc.longitude);
    else
        DEM = false; 
    end
    
    mpolyfile = pickfile('*.mpoly',Inf);
    if opt.cleanup, cellfun(@delete,mpolyfile); mpolyfile = {}; end
    
    cd(basepath);
    for type = opt.types(:)'
        
        if isempty(dir(type{1})), mkdir(type{1}); end
        
        for j = 1:numel(opt.files2copy)
           copyfile(fullfile('./common',opt.files2copy{j}),type{1});
        end
        
        cd(['./' type{1}]);
        
            if opt.cleanup
               delete '*.tpoly'; 
               delete '*.mounts'; 
               delete '*.arrdef';
               delete '*.mpoly';
            end
            samplesystem(type{1},mdims,'landscape',DEM,'files',true,'location',Loc);
            
            lbl{1} = regexprep(pickfile('*.mounts',1),'[\.\\/]*(.*)_[\dxm]+\.mounts','$1');
            lbl{2} = regexprep(pickfile('*.arrdef',1),'[\.\\/]*sample(.*).arrdef','$1');
            prjlbl = [lbl{:} '.ssp'];
                
            if opt.cleanup || isempty(dir(prjlbl))
                copyfile('../common/base.ssp',prjlbl);
                if ~opt.flat, copyfile(fullfile('../common',DEMfile),'./'); end
                f = splitscript(prjlbl,{'phystrck','physmod','layout','terrain','arrdef'},'onerror',opt.onerror);
                assert(~stoponerror || f == 0,'error while evaluatig testbase');
                if ~opt.flat, delete(DEMfile); end
            end
            
            if opt.solve
                f = splitscript(prjlbl,'onerror',opt.onerror);
                assert(~stoponerror || f == 0,'error while evaluatig testbase');
            end

            % Leave only one copy of *.mpoly file, in COMMON folder
            if isempty(mpolyfile)
                mpolyfile = pickfile('*.mpoly',1);
                movefile(mpolyfile,fullfile('../common',mpolyfile));
            else
                cellfun(@delete,pickfile('*.mpoly',Inf));
            end

        cd('..');
    end
end

function F = getlandscape(file,lat,lon,R)
% Return an interpolant following the landscape in geotiff FILE, limited to radius R around
% LAT, LON.

    narginchk(3,4);
    if nargin < 4, R = 200; end
        
    Terr = geotiffread(file);
    
    % Get a WGS84-projected circle of radius R + lattice-size
    c = polygon(60);
    [c.y,c.x] = prj2abs(c.x*R,c.y*R,lat,lon);
    c = offsetpolygon(c,Terr.info.map_info.dx);
    
    % Clip raster to circle limits
    inxbounds = Terr.x > min(c.x) & Terr.x < max(c.x);
    inybounds = Terr.y > min(c.y) & Terr.y < max(c.y);
    
    if ~any(inxbounds) || ~any(inybounds)
       error('Geotiff file does not seem to match the site'); 
    end
    
    [X,Y] = meshgrid(Terr.x(inxbounds),Terr.y(inybounds));
    Z = double(Terr.z(inybounds,inxbounds));
    
    % Project to UTM, offset to lat,lon, and generate interpolant
    [X,Y] = abs2prj(Y,X,lat,lon);
    F = scatteredInterpolant(X(:),Y(:),Z(:),'natural');
end