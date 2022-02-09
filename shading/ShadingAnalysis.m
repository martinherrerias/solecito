function ShRes = ShadingAnalysis(varargin)
% RES = SHADINGANALYSIS(SUNPOS,TRCK,[HOR,SHGEOM,OPT,'opt',val,..])
%   Perform geometrical shading analysis for the tracker field TRCK, at the solar positions
%   given in SUNPOS, optionally considering horizon profile(s) HOR, and using a diffuse shading
%   geometry (SHADINGREGIONS object) SHGEOM.
%
% RES = SHADINGANALYSIS('backup',S) - attempt to resume interrupted calculation from backup
%   structure S = load(FILE), after a previous call with ..,'backup',FILE.
%
% INPUT:
%
%   SUNPOS - Structure with effective azimuth (Az) and elevation (El) degree angles.
% 		Azimuth convention N2E (astronomical).
% 		All fields must be Nt vectors, where Nt is the number of simulation steps.
%
%   TRCK.name - (optional)
%   TRCK.type - {1aF,1aC,1aV,2a,0a} and, depending on TRCK.type: { tilt, azimuth, slope, ..
%       tracklimits, backtracking, and groundcoverratio }. See MOUNTROTATIONS for details.
%   TRCK.origin - [longitude °N,latitude °E, elevation] of project origin
%   TRCK.rotation - Degrees N2W (CCW from North) of project's Y axis.
%   TRCK.geom.border - polygon with tracker outline in Tracker-Coordinate-System (TCS)
%   TRCK.axisoffset - vector from center of rotation to center of tracker, in TCS ([0;0;0])
%   TRCK.centers - 3·Ntr array of center point coordinates, for all TRCK Ntr in the array.
%   TRCK.centerheight - scalar, required for ground-shade calculation
%   TRCK.analysedtrackers - vector of Nu indices for trackers under analysis. (1:Ntr if omitted)
%   TRCK.analysedpoints - 3·Npt array, representing different points over the tracker's surface
% 		(TCS) for which the diffuse shading analysis is to be repeated.
%   TRCK.masks - Nu·Ntr or Ntr·Ntr boolean array, determines which neighboring trackers are 
% 		considered for the shading calculation of each analysed tracker.
%
%   HOR.fHor - polygon3d object for far-horizon-profile (points over the unit sphere)
%   HOR.nHor - Ntr vector of polygon3d objects (near-horizon-profiles for each mount)
%           default is a flat horizon (OPT.flathorizon = true)
% 
%   SHGEOM - SHADINGREGIONS object defining the geometry of diffuse-shading-regions.
%            Default is SHADINGREGIONS() that reads from global SimOptions
% 
%   OPT.sensors - CURRENT: 3xM array of (static) irradiance-sensor normals
%               FUTURE: could include positions, offsets, near-horizons, dynamic rotations, etc.
%
%   OPT.diffuseshading - Otherwise calculate only beam shading    
%   OPT.groundshading - Calculate albedo reduction by projected ground shades (slow)
%   OPT.flathorizon - Ignore horizon profile data
%   OPT.IAM - Incident-Angle-Modifyer function/model (see POLYPROJECTOR and CHECKIAM).
%   OPT.RelTol - relative simulation tolerance
%   OPT.angtol - (degrees) Detail for construction of scenario polygons, recommended sqrt(RelTol)
%   OPT.plotting - Switch between waitbar and graphical display (changed on runtime)
%   OPT.backup - preliminary result file name
%   OPT.autosavemin - dump preliminary results to disk every XX minutes
%
% OUTPUT: SHADINGRESULTS object RES
%
% FUTURE WORK:
%  - Handle multiple mount types / missing modules.
%  - Ground shades over irregular terrain. Projection over triangular mesh, merging and handling of
%    holes in merged 3D shade. Possibly using constrained delaunay triangulation (DelaunayTri).
%  - Beam-shading calculation assumes all mounts have similar orientations (calculation on average
%    plane). A clustering algorithm should be used to allow e.g. building-integrated systems.
%  - For diffuse shading: use a spherical polygon clipping library to merge polygons before
%    projection. Alternatively, consider the possibility of clustering objects into regions
%    with view-field < 90° and using a Gnomonic projection (no edge-refinement required as great-
%    arcs map to straight lines), merging, and reprojecting to IAM.
%
%  - OR... just rewrite everything using a graphics library, without trying to reinvent the wheel.
%
% See also: SHADINGRESULTS, POLYPROJECTOR, MOUNTROTATIONS, CHECKIAM, POAIRRADIANCE

%% Preprocessing

    global SimOptions
    SimOptions = completeoptions(); % fill with defaults if necessary
 
    % Required simulation options
    REQ_OPT = {'flathorizon','diffuseshading','groundshading',...
                'plotting','autosavemin','backup',...
                'RelTol','angtol','cellshadingthreshold','minSunEl','masktolerance'...
                'IAM','sensors'}; % TODO: these two are actually important inputs!
            
    % Mininum list of variables that allow resuming after crash
    BACKUP_VARS = {'sunaz','sunel','Trackers','Terrain','ShRegions','options',...
        'simtimer','BLPoly','BshF','Nsb','DshF','DwF0','DwF','belowhorizon','t'};

    options = rmfield(SimOptions,setdiff(fieldnames(SimOptions),REQ_OPT));
    clear SimOptions
    options.backup = '';
    [options,varargin] = getpairedoptions(varargin,options);

    switch numel(varargin) 
    case 0
    % RES = SHADINGANALYSIS('backup',S)
    % Resume from crash, load existing results from previous partial simulation structure...
        assert(isstruct(options.backup) && all(isfield(options.backup,BACKUP_VARS)),...
            'Failed to recover from backup');

        varargin = struct2cell(orderfields(options.backup,BACKUP_VARS));
        [sunaz,sunel,Trackers,Terrain,ShRegions,options,...
            simtimer,BLPoly,BshF,Nsb,DshF,DwF0,DwF,belowhorizon,t0] = deal(varargin{:});
        clear varargin

        t0 = t0 + 1; % start at the next not-finished timestep
    case {2,3,4,5}
    % RES = SHADINGANALYSIS(SUNPOS,TRCK,[HOR,SHGEOM,OPT])
        if numel(varargin) == 5
            options = completestruct(varargin{end},options);
            options = rmfield(options,setdiff(fieldnames(options),REQ_OPT));
            
            % % Parse options structure: check required, remove unused
            % missing = ~isfield(options,REQ_OPT);
            % assert(~any(missing),'Missing %s',shortliststr(REQ_OPT(missing),'option','quotes',''''));
            % options = cell2struct(cellfun(@(x) options.(x),REQ_OPT','unif',0),REQ_OPT');
        end
        varargin(end+1:4) = {[]};
        [SunPos,Trackers,Terrain,ShRegions] = deal(varargin{1:4});
        clear varargin
        
        % Parse input, PENDING: requires further testing of input

        t0 = 1; % start from the beginning
        simtimer = stopwatch();

        parsestruct(SunPos,{'Az','El'},'-n','-f','-r','-v','-e');
        parsestruct(Trackers,{'analysedpoints','analysedtrackers','axisoffset','centerheight','centers'},'-n','-r','-f');
        assert(isfield(Trackers,'geom') && isa(Trackers.geom,'pvArea'),'Missing/bad TRCK.geom');
        
        sunaz = single(SunPos.Az(:));
        sunel = single(SunPos.El(:));
        clear SunPos
    
        if isempty(ShRegions)
            ShRegions = ShadingRegions();
        end
        assert(isa(ShRegions,'ShadingRegions'),'Invalid SHGEOM');
        
        if ~options.flathorizon && (isempty(Terrain) || ...
                (isempty(Terrain.fHor) && isempty(Terrain.nHor)))
            options.flathorizon = true;
            warning('ShadeAnalysis:flathor','No horizon profile(s)!, a flat horizon will be used');
        end

        if ~options.flathorizon
            switch numel(Terrain.nHor)
                case size(Trackers.centers,2)
                    Terrain.nHor = Terrain.nHor(Trackers.analysedtrackers);
                case numel(Trackers.analysedtrackers) % do nothing
                otherwise
                    error('ShadeAnalysis:nHor','HorProf.nHor must be an Ntr or Nu vector of polygons');
            end
        end
        
        if options.groundshading && ~options.diffuseshading
           warning('Ground-shading requires diffuse-shading! - will be disabled');
           options.groundshading = false;
        end
        
    otherwise
        error('Expecting 0 (crash-resume) or 4-5 positional arguments: SUNPOS,TRCK,HOR,SHGEOM,OPT');
    end
    clear varargin
    
    % Everything works internally in a rotated 'Project Coordinate System', so in order to keep
    % sky patches fixed to an absolute (geographic) reference, we have to work with a rotated
    % version:
    rShReg = rotate(ShRegions,-Trackers.rotation);

    % create an object to handle all plotting functions (now just a waitbar)
    plotobj = shadeanalysisplotter(false,[],sunaz,sunel,[],options);
    lastwill = onCleanup(@() plotobj.close(true)); % close on error
    
    warning_resetter = naptime(); %#ok<NASGU> used to avoid repeated warnings in loop

    % Process 
    Nt = numel(sunaz);
    Ntr = size(Trackers.centers,2);         % all trackers - shadow throwers
    Nu = numel(Trackers.analysedtrackers);  % analysed trackers (receivers)

    Nm = Trackers.geom.dims(1);
    Ne = Trackers.geom.dims(2);
    Nc = Trackers.geom.dims(3);
    
    UINT16MAX = 2^16-1;
    packU16 = @(x) uint16(x*UINT16MAX);
    U16tol = packU16(options.RelTol);
    unpackU16 = @(x) double(x)/UINT16MAX;
    
    % Create a (packed) copy of Trackers.geom, with tighter fits to pv-active-surfaces
    trckgeom = copy(Trackers.geom);
    for j = 1:Nm
        trckgeom.elements(j).border = areaenvelope(trckgeom.elements(j));
    end
    trckgeom.border = areaenvelope(trckgeom);
    trckgeom.border = intersectpolygons(Trackers.geom.border,trckgeom.border);
    trckgeom = struct(trckgeom,'pack16'); % pack16 all border polygons
    trackerarea = polygon.packedarea(trckgeom.border);

    % Flat Horizon (hemishpere) - an approx circle just slightly above XY plane
    FlatHor = polygon3d(linspace(0,360,ceil(1.5*360/options.angtol)),1,eps(pi),'pol');

    % Main projection function is integrated IAM·cos(z)
    IAMprj = polyprojector('IAM','fiam',options.IAM,'normalize',false,...
                                            'angtol',options.angtol,'tol',options.RelTol);
    % IAMprj_r0 = IAMprj.fun(IAMprj.c0);
    
    % Get rotation matrices R, sun vectors, and layout update function
    [TrckFaces,getTrckVertices,sunvec,Rtrck] = plantlayout(Trackers,sunaz,sunel);

    n = permute(Rtrck(:,3,:,:),[1,4,3,2]); % [3,Nt*,Nu*] array of surface normals
    behindsurface = shiftdim(sum(n.*sunvec',1) <= IAMprj.c0,1); % [Nt,Nu*] logical, dot(n,s) <= c0

    equalrotations = size(Rtrck,3) == 1;
    staticsystem = size(Rtrck,4) == 1;
    if staticsystem, Nt_ifmoving = 1; else, Nt_ifmoving = Nt; end
    
    if options.diffuseshading
        % An additional (more restrictive) mask is used for diffuse shading
        mask_diffuse = gettrackermasks(Trackers,options.RelTol,options.IAM);
    
        pointoffsets = Trackers.analysedpoints;
        if size(pointoffsets,1)==2, pointoffsets(3,:) = 0; end
        Np = size(pointoffsets,2);
    else
        Np = 1;
        options.groundshading = false;
    end

    if options.groundshading
        % FUTURE: allow use of external Terrain.nGnd 
        
        % Get ground triangulation(s) and shading update function
        [getgroundshading,Gnd,GndVV,GndP,~] = groundshading(Trackers,sunvec,getTrckVertices);
        
        % make sure most surface normals point up
        Gnd_normals = Gnd.faceNormal;
        k = Gnd_normals(:,3) < 0;
        if mode(k) == 1
            Gnd_normals = -Gnd_normals;
            Gnd.ConnectivityList = Gnd.ConnectivityList(:,[1,3,2]);
            k = ~k;
        end
        if any(k)
            warning('Downward facing surface normals have not been tested');
        end
        
        if staticsystem
            GndVV_W0 = cell(Nu,1);
            GndVV_W = cell(Nu,Np);
        end
    else
        Gnd = [];
    end
    
    if ~isfield(Trackers,'masks')
    % Use as global mask, i.e. also for beam shading
        mask_beam = gettrackermasks(Trackers,options.masktolerance);
    else
        mask_beam = Trackers.masks;
    end

    % Horizon profile(s)
    if options.flathorizon
        Terrain.fHor = FlatHor;
        Terrain.nHor = polygon3d.empty();
    end 

    if t0 == 1
    % Initialize result structures

        % beam shading result will a cell array of polygons (for each timestep and tracker) 
        BLPoly = cell(Nt,Nu);  % int16-packed polygon (see SHADINGRESULTS)

        % See which solar positions lie below the horizon (for each analysed-trackers) 
        if options.flathorizon
            belowhorizon = repmat(sunel < 0,1,Nu);
            Terrain.nmax = 0;
            Terrain.fmax = 0;
        else
        % use the azimuthal projections of both the horizon-profiles and the sun-positions
            AzPrj = polyprojector('azim','clip',false);
            projsp = AzPrj.prj(sunvec');

            AzPrj.angtol = 32/60; % Use sun-diameter of 32" as projection resolution 
            
            projfh = AzPrj.project(Terrain.fHor,[0;0],true);
            for tr = Nu:-1:1
                projnh = AzPrj.project(Terrain.nHor(tr),[0;0],true);
                projnh = mergepolygons(projnh,projfh);
                belowhorizon(:,tr) = ~inpolygon(projsp(1,:),projsp(2,:),projnh.x,projnh.y)';
            end
        end
        
        BshF = zeros(Nt*Nu,Nm,'uint16');
        Nsb = zeros(Nt*Nu,Nm,'uint8');
        BshF(belowhorizon | behindsurface,:) = UINT16MAX; 
        Nsb(belowhorizon | behindsurface,:) = Ne;
        BshF = reshape(BshF,Nt,Nu,Nm);
        Nsb = reshape(Nsb,Nt,Nu,Nm);

        if options.diffuseshading
            % Diffuse-Weight-Factors (Transposition Factors) measure the area of the visible dome 
            % occupied by each region, not considering near-horizons or obstacles.
            [DwF(1:Nu,1).sky] = deal(zeros(Nt_ifmoving,rShReg.n.sky,'uint16'));
            [DwF(1:Nu,1).albedo] = deal(zeros(Nt_ifmoving,rShReg.n.albedo,'uint16'));
            [DwF(1:Nu,1).solar] = deal(zeros(Nt,rShReg.n.solar,'uint16'));

            % Diffuse-Shading-Factors measure the ratio of visible area before and after considering
            % near-horizons and obstacles, for each region.
            [DshF(1:Nu,1:Np).sky] = deal(zeros(Nt_ifmoving,rShReg.n.sky,'uint16'));
            [DshF(1:Nu,1:Np).albedo] = deal(zeros(Nt_ifmoving,rShReg.n.albedo,'uint16'));
            [DshF(1:Nu,1:Np).solar]= deal(zeros(Nt,rShReg.n.solar,'uint16'));
            
            if options.groundshading
                [DwF(1:Nu,1).gndbeam] = deal(zeros(Nt,rShReg.n.albedo,'uint16'));
                [DshF(1:Nu,1:Np).gndbeam] = deal(zeros(Nt,rShReg.n.albedo,'uint16'));
            end
        else
            [F_sky,F_gnd,F_cs] = viewfactors(rShReg,IAMprj,n,sunvec',Terrain.fHor); % [Nt,Nc,Nu*]
            for j = Nu:-1:1   
                DwF(j,1).sky = packU16(F_sky(:,:,j));
                DwF(j,1).albedo = packU16(F_gnd(:,:,j));
                DwF(j,1).solar = packU16(F_cs(:,:,j));
            end
            DshF = []; 
        end
      
        % PROVISIONAL: Get Sensor Transposition Factors, (static & only far-horizon-profile!)
        % TODO: parse a sensor structure that allows tracker-mounted GTI sensors, and static
        %   (shaded) GHI/DHI/GTI sensors at specific locations on the field.
        sensornormals = [0;0;1];
        if ~isempty(options.sensors) 
            
            if istable(options.sensors), options.sensors = table2struct(options.sensors); end
            sensornormals = repmat([0;0;1],1,numel(options.sensors));
            
            for j = 1:numel(options.sensors)
               if isfield(options.sensors,'normal') && ~isempty(options.sensors(j).normal)
                   sensornormals(:,j) = options.sensors(j).normal(:);
               else
                   if isfield(options.sensors,'az'), az = options.sensors(j).az; end
                   if isempty(az), az = 0; end
                   if isfield(options.sensors,'tilt'), el = 90 - options.sensors(j).tilt; end
                   if isempty(el), el = 0; end
                   sensornormals(:,j) = sph2cartV(az,el)';
               end
            end
                            
            % sensornormals = unique(sensornormals','rows')';
            sensornormals = permute(sensornormals,[1,3,2]);
        end
        [F0_sky,F0_gnd,F0_cs] = viewfactors(rShReg,IAMprj,sensornormals,sunvec',Terrain.fHor);  % [Nt,Nc,Ns]
        for j = size(F0_sky,3):-1:1   
            % Sensor diffuse component weight factors. [Ns,1] vector of structures, equivalent to
            % DwF but for Ns sensor orientations instead of Nu mount planes.
            DwF0(j,1).sky = packU16(F0_sky(:,:,j));
            DwF0(j,1).albedo = packU16(F0_gnd(:,:,j));
            DwF0(j,1).solar = packU16(F0_cs(:,:,j));
        end
    end

    AzPrj.clip = 90;

    PjReg0(Nu,1) = struct('horizon',{[]},'sky',{[]},'albedo',{[]},'solar',{[]});
    if staticsystem, saved_pjreg(Nu,Np) = PjReg0(1); end
    trckmask = false(Nu,Ntr);
    partlybehind = false(Nu,Ntr);
    
    simtimer.add2lap('preprocessing');
    % simtimer.resetglobal();
    lastsaved = simtimer.globaltime();
    fibo = [0 1]; % DEBUG: save simtimer splits @ t = 1,2,3,5,8,13...
        
    plotobj.Gnd = Gnd;
    plotobj.Mounts = Trackers;
    % plotobj.Mounts.geom = trckgeom;
    if plotobj.isUI, setbtn(plotobj,[],'Enable','on'); end
    if options.plotting
        plotobj.switchtoplotter(); % switch to full-plotter
        plotobj.pausing = true;
    else
        plotobj.fighandle.reset(); 
    end

%% Time-Step Loop
for t = t0:Nt
%% TIME-LOOP START-UP: Generate scenario (including beam-shade projections)

    if t-t0 == sum(fibo)
        simtimer.addsplit(sprintf('main_loop_%d',sum(fibo)));
        fibo = fibo(2) + [0 fibo(1)];
    end
    if (t == t0+1) && ~options.plotting, plotobj.fighandle.reset(); end
        
    if t ~= t0 && sunel(t) <= options.minSunEl
        continue; 
    end
   
    s = sunvec(t,:)';

    % Indexing variable for Nt(only-if-moving) arrays
    if staticsystem, t_or_1 = 1; else, t_or_1 = t; end

    % Generate scenario: replicate rotated & translated mounts
    if ~staticsystem || t == t0
        if staticsystem, TrckV = getTrckVertices; else, TrckV = getTrckVertices(t); end
        TrckP = polygon3d.vf2poly(TrckV,TrckFaces);
        
        if equalrotations
            Rap = Rtrck(:,:,1,t_or_1);
        else
            % Get an approximate "average-plane rotation" for each timestep
            n = mean(Rtrck(:,3,:,t_or_1),3); % mean surface normal
            Rap = [null(n'),n./norm(n)];

            % TODO: the average plane approach is clearly not very robust when using arrays in very
            % different orientations. Maybe a clustering approach is needed?
            if any(permute(Rtrck(:,3,~behindsurface(t_or_1,:),t_or_1),[3,1,2,4])*Rap(:,3) < 0.5)
                warning('ShadeAnalysis:normals',...
                    'Surface normals in wildly different directions, results might be unreliable');
                warning('off','ShadeAnalysis:normals');
            end
        end
    end

    simtimer.add2lap('scenario');

    % GndShP = polygon.empty; % clear
    if options.groundshading && ~all(belowhorizon(t,:))
        
        [STRI,BTRI,Gnd_shading,Gnd_brightness] = getgroundshading(t);
        Gnd_shading = Gnd_shading.*Gnd_brightness;
    else
        STRI=[]; BTRI = [];
    end
    simtimer.add2lap('ground');

    if ~all(belowhorizon(t,:) | behindsurface(t,:))
    % The projections of all modules in the direction of the sun and over an 'average' tracker 
    % plane are calculated here once and for all (for the current time-step), they will be used
    % later (&) for the beam shading analysis.
        
        % Project in direction s over plane at origin with normal n
        n = Rap(:,3);
        Vp = TrckV - TrckV*n/dot(n,s)*s';
        
        % Change to plane's system of coordinates
        Vp = Vp*Rap;
        
        % Discard z (should be zero), and create packed polygons (for fast clipping)
        POAshP = pack(polygon3d.vf2poly(Vp,TrckFaces));

        simtimer.add2lap('beam');
    end
    if t==t0, simtimer.addsplit('reached_trck_loop'); end

%% Tracker-Loop
    for tr = 1:Nu
    %% TRACKER-LOOP START-UP: get tracker masks and auxiliary coord. system
    
        % GndSh_POA = polygon();
        iatr = Trackers.analysedtrackers(tr); % Index of Analysed TRacker
        
        % Indexing variable for Nu_ifnotequal arrays
        if equalrotations, tr_or_1 = 1; else, tr_or_1 = tr; end

        % Get rotation matrix and offset vector for Point-of-View (POV) coordinate system
        % i.e. inverse of rotation matrix for tracker under analysis
        if equalrotations, Rpov = Rtrck(:,:,1,t_or_1)';
        else, Rpov = Rtrck(:,:,iatr,t_or_1)';
        end

        % get the center of the tracker, including axisoffset
        Pref0 = Trackers.centers(:,iatr)+Rpov'*Trackers.axisoffset(:);

        simtimer.add2lap('beam'); %

        if ~staticsystem || t == t0
        % Get the list of trackers which can actually be 'seen' from the analyzed tracker
            trckmask(tr,:) = mask_beam(tr,:); % start with precalculated mask
            
            if any(trckmask(tr,:))
                for j = find(trckmask(tr,:)) 
                    % Transform every tracker to analysed-tracker coordinate system
                    P = polytranslate(TrckP(j),-Pref0);
                    P = polyrotate(P,Rpov);
                    % check if at least one point lies in front of the tracker plane...
                    trckmask(tr,j) = trckmask(tr,j) && any(P.z > 0); 

                    % avoid projecting beam shades from behind (§)
                    partlybehind(tr,j) = trckmask(tr,j) && any(P.z < 0); 
                end
            end
            simtimer.add2lap('masks');
        end

        %% Beam Shading
        
        if belowhorizon(t,tr) || behindsurface(t,tr_or_1)
            % Sun below horizon or behind tracker-plane
            BshF(t,tr,:) = UINT16MAX; 
            Nsb(t,tr,:) = Ne;
        else
        % Beam Shading Analysis is performed with flat projections over the 'average' plane (&)
        % Any Shading-Tracker projections are substracted from the Analyzed-Tracker's projection,
        % and the final result is then projected back to the analyzed tracker's plane.
        
            if any(trckmask(tr,:))
                % Starting with projection of analysed tracker over 'average' plane...
                % substract any intersecting beam shades (only trackers fully in front, first) (§)
                sunlitPOA = polyclip(POAshP(iatr),POAshP(trckmask(tr,:) & ~partlybehind(tr,:)),0,2,2); % pack64

                if isempty(sunlitPOA)
                % Mount is completely shaded
                    BshF(t,tr,:) = UINT16MAX; 
                    Nsb(t,tr,:) = Ne;
                    
                elseif polygon.packedarea(sunlitPOA) == polygon.packedarea(POAshP(iatr))
                % No shading, BshF(t,tr,:) = Nsb(t,tr,:) = 0
                
                elseif ~isempty(sunlitPOA)
                % Potential shading
                    
                    % Project back into analyzed-tracker's plane
                    P = polyrotate(polygon3d(polygon.unpack(sunlitPOA)),Rap);
                    P = polyprojectshadow(P,s,Rpov(3,:)',Pref0);
                    % ... and go to tracker's system of coordinates
                    P = polytranslate(P,-Pref0);
                    sunlitPOA = pack16(polygon(polyrotate(P,Rpov)));

                    % Clip to PV-active surface
                    sunlitPOA = polyclip(sunlitPOA,trckgeom.border,1); % pack16
                    
                    if ~isempty(sunlitPOA) && any(partlybehind(tr,j))
                    % Substract intersecting beam shades for partially visible trackers (§)
                    
                        for j = find(partlybehind(tr,j))
                            % Clip partially visible polygon above iatr-plane !!!
                            P = clipwithplane(TrckP(j),Rpov(3,:)',Pref0);
                            % Project directly iatr-plane and transform to pov coordinate system
                            P = polyprojectshadow(P,s,Rpov(3,:)',Pref0);
                            P = polytranslate(P,-Pref0);
                            P = polygon(polyrotate(P,Rpov));
                            sunlitPOA = polyclip(sunlitPOA,pack16(P),0,2,2); % pack16
                        end
                    end
                    if ~isempty(sunlitPOA)
                        % PROVISIONAL?: filter out meaningless polygons
                        meaningful = arrayfun(@(p) abs(polygon.packedarea(p))/trackerarea > options.RelTol,sunlitPOA);
                        sunlitPOA(~meaningful) = [];
                    end
                    
                    % Calculate module shading fractions, and number of shaded blocks
                    [BshF(t,tr,:),Nsb(t,tr,:)] = ModuleShadingFactors(sunlitPOA);
                    
                    if ~all(BshF(t,tr,:) == 0) && ~all(BshF(t,tr,:) == UINT16MAX)
                    % 0 < ShF < 1 : partially shaded, save the beam-light-polygon
                        BLPoly{t,tr} = sunlitPOA;
                    end
                end
            else
                % No shading, BshF(t,tr,:) = Nsb(t,tr,:) = 0
            end
        end
        simtimer.add2lap('beam');
        
        plotobj.refresh(); % apply any queued request to change plotter instance
        if plotobj.plotting && (options.diffuseshading || isstruct(BLPoly{t,tr}))
            
            % Plotting of Scenario & Beam Shading
            if isempty(BLPoly{t,tr}) && all(BshF(t,tr,:)==0)
                plotobj.majorupdate(t,tr,1,STRI,BTRI);
            else
                plotobj.majorupdate(t,tr,BLPoly{t,tr},STRI,BTRI);
            end

            if ~options.diffuseshading
                drawnow(); % last thing to plot
                simtimer.add2lap('plotting');
                if plotobj.pausing, plotobj.pause(); simtimer.resetlocal(); end
            end
            simtimer.add2lap('plotting');
        end
        if ~options.diffuseshading
            if ~plotobj.plotting, plotobj.updatewaitbar((t-t0)/(Nt-t0+1),'-addtime'); end
            continue; % skip the rest
        end

        %% Region-Projection & Transposition Factors (far-horizon only)
        if (~equalrotations || tr == 1)
            if (~staticsystem || t == t0)
            % For each mount i.e. surface orientation...
                
                % Calculate projections of static regions
                [prj0,~,prj0_NC] = projectstatic(rShReg,IAMprj,Terrain.fHor,Rpov');

                % Get [possibly provisional (#)] sky view factors
                DwF(tr).sky(t_or_1,:) = packU16(cellfun(@area,prj0.sky)./rShReg.solidangles.sky);

                % Pack16 (significant) polygons, for faster clipping later
                visible_sky = DwF(tr).sky(t_or_1,:) >= U16tol;
                DwF(tr).sky(t_or_1,~visible_sky) = 0;
                PjReg0(tr).sky = cell(1,rShReg.n.sky);
                PjReg0(tr).sky(visible_sky) = cellfun(@pack16,prj0.sky(visible_sky),'unif',0);
                % PjReg0(tr).sky(~visible_sky) = {struct.empty};
                
                % Get ground view factors and pack16 (common) horizon projection...
                PjReg0(tr).horizon = pack16(prj0.horizon);

                DwF(tr).albedo(t_or_1,:) = packU16(cellfun(@area,prj0.albedo)./rShReg.solidangles.albedo);

                visible_gnd = DwF(tr).albedo(t_or_1,:) >= U16tol;
                DwF(tr).albedo(t_or_1,~visible_gnd) = 0;
                
                if options.flathorizon
                    PjReg0(tr).albedo = cell(1,rShReg.n.albedo);
                    PjReg0(tr).albedo(visible_gnd) = cellfun(@pack16,prj0.albedo(visible_gnd),'unif',0);
                    % PjReg0(tr).albedo(~visible_gnd) = {struct.empty};
                else
                % ... keep non-clipped prj0_NC.albedo to use (higher) near-horizons later (#)
                    PjReg0(tr).albedo = cellfun(@pack16,prj0_NC.albedo,'unif',0);
                end
            else
                % Use previously calculated PjReg0(tr), solar regions will be updated anyway
            end
            
            % Calculate projections & view-factors for circumsolar regions
            PjReg0(tr).solar = projectsolar(rShReg,IAMprj,s,Rpov')';
            PjReg0(tr).solar = cellfun(@pack16,PjReg0(tr).solar,'unif',0);
            % PjReg0(tr).solar = cellfun(@(p) polyclip(p,PjReg0(tr).horizon,'i'),PjReg0(tr).solar,'unif',0);
            
            cs_areas = cellfun(@polygon.packedarea,PjReg0(tr).solar);
            cs_areas(:,2:end) = diff(cs_areas,1); % ring areas
            DwF(tr).solar(t,:) = packU16(cs_areas./rShReg.solidangles.solar);    
            visible_cs = DwF(tr).solar(t,:) >= U16tol;
            DwF(tr).solar(t,~visible_cs) = 0;
            PjReg0(tr).solar(~visible_cs) = {struct.empty};
        else
            PjReg0(tr) = PjReg0(1);
            DwF(tr) = DwF(1);
        end
        
        % if ~options.diffuseshading
        %     plotter.refresh(); % apply any queued request to change plotter instance
        %     if ~plotter.plotting && mod(tr-1,100) == 0
        %         plotter.updatewaitbar((tr+(t-t0)*Nu)/((Nt-t0+1)*Nu),'-addtime');
        %     end
        %     continue; % SKIP TO NEXT TRACKER: EVERYTHING BELOW ONLY IF DIFFUSESHADING!
        % end
        % simtimer.add2lap('plotting');
        
        %% Diffuse Shading: near horizons
        if ~options.flathorizon

            if (~staticsystem || t == t0)
            % Apply Local Horizon profiles (#)
                P = IAMprj.project(polyrotate(Terrain.nHor(tr),Rpov),Rpov(3,:)',true);
                PjReg0(tr).horizon = polyclip(PjReg0(tr).horizon,pack16(P),'i');
                
                visible_sky = ~cellfun(@isempty,PjReg0(tr).sky);
                PjReg0(tr).sky(visible_sky) = cellfun(@(p) ...
                    polyclip(p,PjReg0(tr).horizon,'i'),PjReg0(tr).sky(visible_sky),'unif',0);
                DwF_sky_near = packU16(...
                    cellfun(@polygon.packedarea,PjReg0(tr).sky)./ShRegions.solidangles.sky);
                visible_sky = DwF_sky_near >= U16tol;   
                PjReg0(tr).sky(~visible_sky) = {struct.empty};
                % DwF_sky_near(~visible_sky) = 0;
                % DwF(tr).sky(t,:) = DwF_sky_near;     
                for j = 1:Np, DshF(tr,j).sky(t_or_1,:) = DwF_sky_near; end
                
                visible_gnd = ~cellfun(@isempty,PjReg0(tr).albedo);
                PjReg0(tr).albedo(visible_gnd) = cellfun(@(p) ...
                    polyclip(p,PjReg0(tr).horizon,'d'),PjReg0(tr).albedo(visible_gnd),'unif',0);
                DwF_gnd_near = packU16(...
                    cellfun(@polygon.packedarea,PjReg0(tr).albedo)./ShRegions.solidangles.albedo);
                visible_gnd = DwF_gnd_near >= U16tol;
                PjReg0(tr).albedo(~visible_gnd) = {struct.empty};
                % DwF_gnd_near(~visible_gnd) = 0;
                % DwF(tr).albedo(t,:) = DwF_gnd_near;
                for j = 1:Np, DshF(tr,j).albedo(t_or_1,:) = DwF_gnd_near; end
            end
            PjReg0(tr).solar = cellfun(@(p) ...
                polyclip(p,PjReg0(tr).horizon,'i'),PjReg0(tr).solar,'unif',0);
            
            cs_areas = cellfun(@polygon.packedarea,PjReg0(tr).solar);
            cs_areas(:,2:end) = diff(cs_areas,1); % ring areas
            DwF_cs_near = packU16(cs_areas./rShReg.solidangles.solar); 
            
            visible_cs = DwF_cs_near >= U16tol;
            PjReg0(tr).solar(~visible_cs) = {struct.empty};
            % DwF(tr).solar(t,~visible_cs) = 0;
            % DwF(tr).solar(t,:) = DwF_cs_near;
            for j = 1:Np, DshF(tr,j).solar(t_or_1,:) = DwF_cs_near; end
        end
        simtimer.add2lap('diffuse');
            
        % Determine visible ground shades
        if options.groundshading && ~belowhorizon(t,tr) && any(visible_gnd)
            % GndSh_POA = polytranslate(GndShP,-Pref0);
            % GndSh_POA = polyrotate(GndSh_POA,Rpov);
            % GndSh_POA = clipwithplane(GndSh_POA,[0;0;1],[0;0;0]);
            
            if ~staticsystem || isempty(GndVV_W0{tr})
                W = GndVV_weights(Gnd,GndVV,GndP,Gnd_normals,Pref0,Rpov,IAMprj,PjReg0(tr).albedo);
                GndVV_W0{tr} = W;
            else
                W = GndVV_W0{tr};
            end

            DwF_albedo = unpackU16(DwF(tr).albedo(t_or_1,:));
            DwF(tr).gndbeam(t,:) = packU16((Gnd_brightness'*W).*DwF_albedo);
            
            % assume any areas not "assigned" are flat unshaded ground
            f = ~any(W,1) & visible_gnd;
            DwF(tr).gndbeam(t,f) = packU16(max(0,s(3).*DwF_albedo(f)));
        end
        simtimer.add2lap('ground');

        % Plot POV-independent projected regions
        if plotobj.plotting
            plotobj.plotskyregions(PjReg0(tr),false,...
                DwF(tr).sky(t_or_1,:),DwF(tr).albedo(t_or_1,:),DwF(tr).solar(t,:));
            simtimer.add2lap('plotting');
        end

        if t==t0 && tr == 1, simtimer.addsplit('reached_pt_loop'); end

%         %% Consider near shades (POV-dependent)
        for pt = 1:Np    
            % Move the reference point over the tracker, if required
            Pref = Pref0 + Rpov'*(pointoffsets(:,pt));

            % Map visible trackers to radiance preserving projection
            if (~staticsystem || t == t0) || plotobj.plotting
                
                % Local projection (POV-dependent), start with near-horizon-shaded PjReg0(tr)
                pjreg = PjReg0(tr); 
            
                relevant = trckmask(tr,:) & mask_diffuse(tr,:);
                P = polytranslate(TrckP(relevant),-Pref);         % Go to POV coordinate system
                P = polyrotate(P,Rpov);                           % ...                 
                PjTrck = arrayfun(@IAMprj.project,P,'unif',0);    % get flat projections
                if ~isempty(PjTrck)
                    PjObstacles = pack16([PjTrck{:}]);                          % ... pack16
                    PjObstacles = polyclip(PjObstacles,struct.empty,'u','pos'); % ... and merge
                else
                    PjObstacles = struct.empty();
                end

                if ~isempty(PjObstacles)
                    % Substract the projection of all obstacles from (visible) static sky-regions
                    visible_sky = DwF(tr).sky(t_or_1,:) > 0;
                    pjreg.sky(visible_sky) = cellfun(@(p) polyclip(p,PjObstacles,'d'),pjreg.sky(visible_sky),'unif',0);
                    pjreg.sky(~visible_sky) = {struct.empty};
                    
                    f = cellfun(@polygon.packedarea,pjreg.sky(visible_sky))./rShReg.solidangles.sky(visible_sky);
                    DshF(tr,pt).sky(t_or_1,visible_sky) = packU16(f);
                    
                    visible_sky = visible_sky & DshF(tr,pt).sky(t_or_1,:) > 0;
                    pjreg.sky(~visible_sky) = {struct.empty};
                    
                    % ... and from (visible) albedo-regions
                    visible_gnd = DwF(tr).albedo(t_or_1,:) > 0;
                    pjreg.albedo(visible_gnd) = cellfun(@(p) polyclip(p,PjObstacles,'d'),pjreg.albedo(visible_gnd),'unif',0);
                    f = cellfun(@polygon.packedarea,pjreg.albedo(visible_gnd))./rShReg.solidangles.albedo(visible_gnd);
                    DshF(tr,pt).albedo(t_or_1,visible_gnd) = packU16(f);
                    
                    visible_gnd = visible_gnd & DshF(tr,pt).albedo(t_or_1,:) > 0;
                    pjreg.albedo(~visible_gnd) = {struct.empty};
                else
                    DshF(tr,pt).sky(t_or_1,:) = DwF(tr).sky(t_or_1,:);
                    DshF(tr,pt).albedo(t_or_1,:) = DwF(tr).albedo(t_or_1,:);
                end
                
                if staticsystem
                    saved_pjreg(tr,pt) = pjreg;
                end
            else
            % Load the projection calculated on first run
                pjreg = saved_pjreg(tr,pt);
                pjreg.solar = PjReg0(tr).solar;
            end
            
            % Substract the projection of all obstacles from circumsolar-regions
            visible_cs = DwF(tr).solar(t,:) >= U16tol;
            if any(visible_cs)
                pjreg.solar(visible_cs) = cellfun(@(p) polyclip(p,PjObstacles,'d'),pjreg.solar(visible_cs),'unif',0);
                cs_areas = cellfun(@polygon.packedarea,pjreg.solar);
                cs_areas(:,2:end) = diff(cs_areas,1); % ring areas
                DshF(tr,pt).solar(t,:) = packU16(cs_areas./rShReg.solidangles.solar); 
                visible_cs = visible_cs & DshF(tr,pt).solar(t,visible_cs) >= U16tol;
                pjreg.solar(~visible_cs) = {struct.empty};
            end

            if options.groundshading && ~belowhorizon(t,tr) && any(visible_gnd)

                if ~staticsystem || isempty(GndVV_W{tr,pt})
                    W = GndVV_weights(Gnd,GndVV,GndP,Gnd_normals,Pref,Rpov,IAMprj,pjreg.albedo);
                    GndVV_W{tr,pt} = W;
                else
                    W = GndVV_W{tr,pt};
                end
                
                DwF_albedo = unpackU16(DshF(tr,pt).albedo(t_or_1,:));
                DshF(tr,pt).gndbeam(t,:) = packU16((Gnd_shading'*W).*DwF_albedo);
                
                f = ~any(W,1) & DwF(tr).gndbeam(t,:) > 0;
                DshF(tr,pt).gndbeam(t,f) = DwF(tr).gndbeam(t,f);
                
            end

            % Plot projections
            if plotobj.plotting
                simtimer.add2lap('diffuse');
                
                % plotobj.plotaxis(Pref,Rpov','r-',s); % Plot POV axis
                plotobj.minorupdate(t,tr,pointoffsets(:,pt));

                % Plot (or update) projections of trackers [and ground shades]
                %plotter.plotprojectedtrck(PjTrckP,PjGndShP(shmask));
                % plotobj.plotprojectedtrck(PjTrck,PjShadesGnd);
                
                shf = DshF(tr,pt).albedo(t_or_1,:);
                if options.groundshading
                   shf = double(shf)./UINT16MAX.*(0.2 + 0.8*double(DshF(tr,pt).gndbeam(t,:))/max(0.05,s(3)));
                end

                % Plot shaded sky regions
                plotobj.plotskyregions(pjreg,true,...
                    DshF(tr,pt).sky(t_or_1,:),shf,DshF(tr,pt).solar(t,:)); 
                %plotter.plotskyregions(PVprj_sky,PVprj_sun,PVprj_gnd);
                
                title(sprintf('pt = %d/%d',pt,Np));
                drawnow; 
                    % plotter.exportfigures(t,tr,pt);
                simtimer.add2lap('plotting');

                if plotobj.pausing % pause if required, stop the timer if so
                    plotobj.pause();
                    simtimer.resetlocal();
                end
            end
        end % point-of-view loop
        
        if ~plotobj.plotting && (t == t0 || plotobj.worthupdating)
            if t == t0
                msg = sprintf('Static factors: mount %d/%d...',tr,Nu);
                plotobj.updatewaitbar(tr/Nu,msg,'-addtime');
            else
                msg = sprintf('Dynamic factors: t = %d/%d, mount %d/%d...',t-t0,Nt-t0+1,tr,Nu);
                plotobj.updatewaitbar((tr+(t-t0)*Nu)/((Nt-t0+1)*Nu),msg,'-addtime');
            end
        end
    end % tracker loop
    simtimer.add2lap('diffuse');
                
    % Dump complete workspace to partial-result file every AUTOSAVEMIN minutes
    if ~isempty(options.backup) && ((t == Nt) || ...
            (simtimer.globaltime - lastsaved > options.autosavemin*60))
        save(options.backup,BACKUP_VARS{:}); 
        lastsaved = simtimer.globaltime;
        simtimer.add2lap('backup');
    end
end % time-step loop

%if ~plotter.plotting,(delete(plotter.fighandle)); end 
delete(lastwill); % plotter.close(true);

% Reduce Beam Shading Polygons to a list of not-emtpy elements
someshading = cellfun(@isstruct,BLPoly);
fullshading = all(BshF == UINT16MAX,3);

if ~isempty(options.backup)
    save(options.backup,BACKUP_VARS{:}); % DEBUG
    simtimer.add2lap('backup');
end

%% Pack results to go
ShRes.info = struct('platform',system_dependent('getos'),'code',getGitInfo());
ShRes.info.timing = simtimer;
ShRes.info.options = options;
ShRes.worldgeom = ShRegions;
ShRes.mountgeom = trckgeom;
ShRes.pts = Trackers.analysedpoints;
ShRes.az = sunaz;
ShRes.el = sunel;
ShRes.partshaded = someshading;
ShRes.fullshaded = fullshading;

ShRes.BLPoly = BLPoly(someshading(:));

ShRes.BshF = BshF;
ShRes.Nsb = Nsb;
ShRes.DwF = DwF;
ShRes.DshF = DshF;
ShRes.DwF0 = DwF0;

ShRes = ShadingResults(ShRes);

simtimer.add2lap('packing');
%%
	
function [Fgs,Nsb] = ModuleShadingFactors(LP)
% Calculate linear geometrical shading factors FGS and number of shaded blocks NSB for each module
% of a given the tracker (trckgeom) given light-polygon LP
    
    THRESH = UINT16MAX/Nc*options.cellshadingthreshold;

    Fgs = zeros(1,Nm,'uint16');
    Nsb = zeros(1,Nm,'uint8');
    
    if isempty(LP), Fgs(:) = UINT16MAX; Nsb(:) = Ne; return; end  % No light (ShF = 1, Nsb = Ne)
    if isequal(LP,trckgeom.border), return; end % No shade (Fgs = Nsb = 0)
            
    [sf,LP] = shadingfactor(LP,trckgeom.border,trckgeom.area);
    if sf < THRESH/(Nm*Ne)                                       % No shade (Fgs = Nsb = 0)
    elseif UINT16MAX - sf < THRESH/(Nm*Ne)                       % No light
        Fgs(:) = UINT16MAX; 
        Nsb(:) = Ne;  
    else                                                            % Partial shading?
        for i = 1:Nm
            [Fgs(i),p] = shadingfactor(LP,trckgeom.elements(i).border,trckgeom.elements(i).area);
            if Fgs(i) < THRESH/Ne
            % No shading
            	Fgs(i) = 0;
                Nsb(i) = 0;
            elseif UINT16MAX - Fgs(i) < THRESH/Ne
            % Full shading
                Fgs(i) = UINT16MAX;
                Nsb(i) = Ne;
            else
            % Partial shading
                Nsb(i) = sum(cellfun(@(g,a) shadingfactor(p,g,a) > THRESH,...
                    {trckgeom.elements(i).elements.border},{trckgeom.elements(i).elements.area}));
            end
        end
    end    
    
    function [f,LP] = shadingfactor(LP,G,a)
    % Fraction WITHOUT light-polygon LP on pvArea: f = [1-area(G & P)/area(G)]·MAX_F
        
        LP = polyclip(G,LP,1); % G.border & LP
        f = packU16(1 - polygon.packedarea(LP)/a);
    end
end

end
