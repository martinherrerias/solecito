classdef ShadingRegions
% Object class to represent a unit-sphere polygonal tesellation meant for arbitrary irradiance 
% transposition and diffuse-shading calculations. Each object has three main properties, each a 
% cell-array of spherical polygons (POLYGON3D objects):
%
%   OBJ.sky - static regions that define regions of visible sky (to be clipped above horizon)
%   OBJ.albedo - static regions that define regions of visible ground (to be clipped below horizon)
%   OBJ.solar - moving (i.e. sun-centered) circumsolar regions, in sun-centered coordinates*.
%
% Since each set of regions (sky/albedo) is meant to cover the respective half-dome, and to avoid
% redundancy in calculations, there can be one 'implicit' region for each of these fields, which 
% is simply defined as any area not included in the explicit spherical polygons. As an example,
% the Perez 1990 geometrical framework could be defined as follows:
%
%   OBJ.solar = {polygon3d(0:5:360,90 - 25,1,'sph')}; - 25° circumsolar circle 
%   OBJ.sky = {polygon3d(0:360,10,1,'sph')}; - single explicit region (isotropic sky)
%   OBJ.albedo = {}; - no explicit albedo
%
%   OBJ.implicit.sky - returns the missing, implicit sky region (10° horizon brightening band)
%   OBJ.implicit.albedo - returns the complete ground below horizon
%
% Overall, there are four transient properties beside the three main explicit polygonal regions:
%
%   OBJ.hasimplicit - {.sky,.solar,.albedo} - scalar boolean: do regions cover everything?
%   OBJ.implicit - {.sky,.solar,.albedo} - 3D polygon for areas not already included
%   OBJ.n - {.sky,.solar,.albedo} - number of regions (including implicit)
%   OBJ.haszenith - {.sky,.solar,.albedo} - OBJ.n boolean vectors: is point [0,0,1] in projection?
%
% The last is meant to be used by POLYPROJECTOR.PROJECT to avoid negative projections.
%
% TODO: spherical polygons should be stored using a face-vertex representation, since most
%   vertices are shared between at least two regions.
%
% See also: SHADINGREGIONS.SHADINGREGIONS, PROJECTSTATIC, PROJECTSOLAR, SHADINGREGIONS.PLOT

    
properties (SetAccess = protected, GetAccess = public)
    info    % free structure field (usually with info.name = char key)
    sky     % cell-array of 3D polygons (static sky regions - must include lowest horizon!)
    albedo  % cell-array of 3D polygons (static ground regions - must include highest horizon!)
    cs      % [sorted] vector of circumsolar-region angles (degrees)
end
properties (Transient)
   solar        % cell-array of 3D polygons (moving circumsolar regions - sun-centered coords.)
   hasimplicit  % {.sky,.solar,.albedo} - scalar boolean values: do regions cover everything?
   haszenith    % {.sky,.solar,.albedo} - boolean values for each 3D polygon: [0,0] in projection?
   implicit     % {.sky,.solar,.albedo} - 3D polygon for areas not already included
   n            % {.sky,.solar,.albedo} - number of regions (including implicit)
   solidangles         
end
properties (Constant)
   TINYARC = 0.1 ;   % (degrees) ~1 mm precision for ± 32 m mount-vertex coordinates
end
methods
    function obj = ShadingRegions(varargin)
    % OBJ = SHADINGREGIONS(...) - Creates a unit-sphere tesellation (SHADINGREGIONS object) from
    %   either a predefined named keyword, or explicit cell-arrays of 3D-polygons.
    %
    % OBJ = SHADINGREGIONS() - Return an empty object (nested empty fields). Note that you must
    %   then use SHADINGREGIONS.LOADOBJ(OBJ) to get a working (isotropic) object.
    %
    % OBJ = SHADINGREGIONS('auto'/'default') - set X = SimOptions.skyregions, and...
    %   a) If X is a KEY or cell-array of 3d-polygons, return SHADINGREGIONS(X)
    %   b) If X is already a SHADINGREGIONS object, return X.
    %   c) If X is a cell array of arguments, name-value pairs, etc. use SHADINGREGIONS(X{:})
    %
    % OBJ = SHADINGREGIONS(KEY) - use one of a series of predefined tesellations:
    %
    %   'isotropic' - Isotropic sky and albedo, no circumsolar regions.
    %   'haydavies' - CS* radius circumsolar region, isotropic sky and albedo.
    %   'perez' - 90° - HBB* isotropic sky (implicit horizon-brightening) + 'haydavies'
    %   'perezN' - Each Perez region is divided into N cardinal directions. Unlike 'perez', it 
    %       also allows vector-valued HBB and/or CS radius (see below).
    %   'sunto' - 11 sky regions reflecting Sunto sensor geometry.
    %   'sunto2' - 21 sky regions reflecting Sunto sensor geometry.
    %   'bucky' - 71 quasi-uniform sky regions with icosahedral-dodecahedral symmetry.
    %   'eko' - 145 zones of the EKO MS-321LR Sky Scanner
    %   'dumortier' - 13 zones of Dumortier (1995)
    %   'tinN' - for N = 40,160,640,2560,..(10·4^z with z > 0) Recursive division of the faces 
    %       of an icosahedral tessellation.
    %   'rbN' - Recursive Bisection N, starting with N approx. square elements on the horizon, and
    %       doubling their size in successive steps towards the zenith. It's recommended that N 
    %       is a multiple of a large enough power of 2 (8,16,32,..), as division will stop when 
    %       bisection is no longer possible.
    %       Albedo is split into N/2 wedge-shaped regions.
    %   'urN' - Quasi-Uniform Random N, using N voronoi cells centered at N points distributed
    %       evenly over the hemisphere (see SPHEREPOINTS).
    %       Albedo is split into N/2 wedge-shaped regions.
    %
    % OBJ = SHADINGREGIONS(S,A) - provide separate arguments for Sky and Albedo regions.
    %   Each argument can be a KEY (e.g. to combine features from two named models), or
    %   a cell-array of explicit 3D-polygons, e.g. the output of MULTISENSOR for divisions
    %   following a spherical-voronoi tesellation.
    %
    %   If omitted (or empty), the resulting SHADINGREGIONS object will have a single (implicit)
    %   region that covers the respective area.
    %
    % OBJ = SHADINGREGIONS(..,'CS',R) - Circumsolar region radii (degrees), passed to OBJ.CS.
    %   R can be a vector, from which concentric rings will be created for OBJ.SOLAR regions.
    %   If omitted (or empty), a single (implicit) 32' diameter sun disc is considered.
    %   (*) For Perez and Hay-Davies models, CS defaults to SimOptions.CSradius.
    %   
    % OBJ = SHADINGREGIONS(..,'HBB',H) - (scalar, for 'perez', vector for 'perezN') horizon-
    %   brightening-band angle(s), in degrees. (*) HBB defaults to SimOptions.HBBwidth.
    %
    % OBJ = SHADINGREGIONS(..,'angtol',TOL) - controls polygon resolution
    % OBJ = SHADINGREGIONS(..,'anisotropic',false) - ignore all arguments, return isotropic object.
    % OBJ = SHADINGREGIONS(..,'diffusemodel',key) - override the meaning of key: 'auto'
    %
    % See also: SHADINGREGIONS, SHADINGREGIONS.PLOT, MULTISENSOR
    
        if nargin == 0
            foo = struct('sky',[],'solar',[],'albedo',[]);
            obj.info = struct(); 
            obj.sky = {};
            obj.albedo = {};
            obj.cs = [];
            obj.hasimplicit = foo;
            obj.haszenith = foo;
            obj.implicit = foo;
            obj.n = foo;
            obj.solidangles = foo;
           return;
        end
    
        ispolygoncell = @(x) iscell(x) && all(cellfun(@(p) isa(p,'polygon3d'),x));

        opt = {'angtol','anisotropic','skyregions','HBBwidth','CSradius'}';
        opt = cell2struct(cellfun(@getSimOption,opt,'unif',0),opt);
        opt.HBB = NaN; opt.CS = NaN;  % HBB, CS are actual user options, 
                                      % HBBwidth, CSradius are (hidden) defaults for Perez
        
        if nargin == 1 && ischar(varargin{1}) && any(strcmpi(varargin{1},{'auto','default'}))
            if isa(opt.skyregions,'ShadingRegions')
            % If SimOptions.skyregions is already an object, make A = ShadingRegions(A)
                obj = opt.skyregions;
            elseif iscell(opt.skyregions) && ~ispolygoncell(opt.skyregions)
                obj = ShadingRegions(opt.skyregions{:});
            else
                assert(~isempty(opt.skyregions) && ...
                    ~(ischar(opt.skyregions) && any(strcmpi(opt.skyregions,{'auto','default'}))),...
                    'Invalid SimOptions.skyregions: causes infinite recursion')
                obj = ShadingRegions(opt.skyregions);
            end
            return;
        end
                           
        [opt,args] = getpairedoptions(varargin,opt);
        if isempty(args), args = {'iso'}; end
        
        if ~opt.anisotropic
            for k = 1:numel(args)
                if (ischar(args{k}) && ~any(strcmpi(args{k},{'auto','default','iso','isotropic'}))) || ...
                   (ispolygoncell(args{k}) && ~isempty(args{k}))
                   warning('Diffuse shading is set to isotropic: ignoring Shading-Region arguments') 
                end
                
                obj = ShadingRegions();
                obj.info.name = 'Isotropic';
                obj = ShadingRegions.loadobj(obj);
                return;
            end
        end

        if numel(args) == 1 && ischar(args{1})  
        % SHADINGREGIONS(key)
            [obj,opt] = namedregions(args{1},opt);
        else
        % SHADINGREGIONS(S,A)
            assert(all(cellfun(@(x) ischar(x) || ispolygoncell(x),args)),'Unrecognized argument(s)');
            %optargs = [fieldnames(opt),struct2cell(opt)]';

            if ischar(args{1})
            % Build from keyword(s)
                [obj,opt] = namedregions(args{1},opt);
            elseif ispolygoncell(args{1})
            % Build from explicit polygon3d tessellation
                obj.info.name = 'explicit';
                obj.sky = args{1};
                obj.albedo = {};
            end
            
            if numel(args) > 1
            % Edit with albedo tessellation, if available
                if ischar(args{2})
                    obj2 = namedregions(args{2},opt);
                    obj.albedo = obj2.albedo;
                elseif ispolygoncell(args{1})
                    obj2.info = 'explicit';
                    obj.albedo = args{2};
                end
                obj.info.name = sprintf('sky: %s / albedo: %s',obj.info.name,obj2.info.name);
            end
        end
        
        % Make sure no region contains zenith or nadir as points
        AzPrj = polyprojector('azim','clip',180-ShadingRegions.TINYARC,'angtol',120);
        [~,obj.sky] = cellfun(@AzPrj.project,obj.sky,'unif',0);
        [~,obj.albedo] = cellfun(@AzPrj.project,obj.albedo,'unif',0);
        AzPrj.south = true;
        [~,obj.sky] = cellfun(@AzPrj.project,obj.sky,'unif',0);
        [~,obj.albedo] = cellfun(@AzPrj.project,obj.albedo,'unif',0);
                
        if isnan(opt.CS), opt.CS = []; end
        obj.cs = sort(opt.CS(:))';
        assert(all(obj.cs > 32/60 & obj.cs < 90),'Invalid circumsolar-angles');
        
        obj = ShadingRegions.loadobj(obj);
        

        function [obj,opt] = namedregions(arg,opt)
        % Return predefined tessellations from keyword ARG
            
            % Start with empty regions
            obj.info.name = arg;
            obj.info.angtol = opt.angtol;
            obj.sky = {}; 
            obj.albedo = {};
            
            % Split keys of the form 'XXNN' into key = 'XX', arg = NN
            errmsg = ['Unknown key: ',arg];
            key = regexp(lower(arg),'^(\D+)(\d*)$','tokens');
            assert(~isempty(key),errmsg);
            [key,N] = deal(key{1}{:});
            if ~isempty(N), N = str2double(N); end
            
            da = @(z) 360/ceil(360/min(0.8*opt.angtol/sind(z),45));
            cregions = @(zz) arrayfun(@(z) polygon3d(da(z)/2:da(z):360,90-z,1,'sph'),zz,'unif',0);

            switch key
            case {'iso','isotropic'}
                assert(isempty(N),errmsg);
                obj.info.name = 'isotropic';
            case 'haydavies'
                assert(isempty(N),errmsg);
                % obj.sky = {} i.e. implicit single static sky region (no horizon brightening)
                % obj.albedo = {} i.e. implicit isotropic ground region
                if isnan(opt.CS), opt.CS = opt.CSradius; end
                obj.info.name = sprintf('Hay-Davies (%0.1f° CSR)',opt.CS);
                obj.info.CS = opt.CS;
            case 'perez'
                % Get Perez's defaults from SimOptions
                if isnan(opt.CS), opt.CS = opt.CSradius; end
                if isnan(opt.HBB), opt.HBB = opt.HBBwidth; end
                isnice = @(x) isscalar(x) && isfinite(x) && x > 0;
                obj.info.CS = opt.CS;
                obj.info.HBB = opt.HBB;
                
                if isempty(N)
                    % Circumsolar region, typically 25°-half angle
                    assert(isnice(opt.CS),'Perez model requires scalar CS radius >= 0');                        
                    % Limit of isotropic sky: upper edge of the horizon brightening band 6.5°
                    assert(isnice(opt.HBB),'Perez model requires scalar > 0 HBB width'); 
                    obj.sky = cregions(90-opt.HBB);
                    obj.info.name = sprintf('Perez et al. (%0.1f° CSR, %0.1f° HBB)',opt.CS,opt.HBB);
                    
                    % obj.albedo = {} i.e. implicit isotropic ground region
                else
                % perezN: split isotropic, HB and albedo into N wedges
                    [~,wedges] = multisensor(90,N/2:360/N:360);
                    obj.albedo = wedges;
                    obj.sky = {};
                    for hbw = opt.HBB
                        HB = [0;0;sind(hbw)];
                        obj.sky = [obj.sky, cellfun(@(p) clipwithplane(p,[0;0;-1],HB),wedges,'unif',0)];
                        wedges = cellfun(@(p) clipwithplane(p,[0;0;1],HB),wedges,'unif',0);
                    end
                    obj.sky = [obj.sky, wedges];
                    obj.info.N = N;
                    obj.info.name = sprintf('perez%d',N);
                end
            case {'eko','MS-321LR'}
            % 145 zones of MS-321LR Sky Scanner
            % https://media.eko-eu.com/assets/media/MS-321LR_Manual.pdf
                obj.info.name = 'EKO MS-321LR Sky-Scanner';
                [~,obj.sky] = multisensor(...
                            90-6,90:-12:-270+12,...  % From azimuth 0° Ch1 to Ch30, clockwise 
                            90-18,90+12:12:450,...   % From azimuth -12° Ch31 to Ch60, counterclockwise
                            90-30,90:-15:-270+15,... % From azimuth 0° Ch61 to Ch84, clockwise
                            90-42,90+15:15:450,...
                            90-54,90:-20:-270+20,... % ...
                            90-66,90+30:30:450,...
                            90-78,90:-60:-270+60,... % From azimuth 0° Ch139 to Ch144, clockwise
                            0,0);                    % Zenith Ch145
                
            case {'satellight','satel-light','dumortier'}
            % 13 zones of Dumortier (1995): http://www.satellight.com/indexgA.htm
                obj.info.name = 'dumortier';
                obj.sky = rings([0 5/13 sqrt(105)/13 1],[1 4 8],opt.angtol*pi/180);
                obj.sky = obj.sky([1,2,5,4,3,7,6,13:-1:8]);     % reorder to match convention
                [~,obj.albedo] = multisensor(90,67.5:-45:-270);
            case 'sunto'
                if isempty(N)
                % sunto: 11 sky regions reflecting Sunto sensor geometry, 5 albedo wedges
                    obj.info.name = 'sunto';
                    [~,obj.sky] = multisensor(0,0,45,(0:4)*360/5,90,(0.5:4.5)*360/5);
                    [~,obj.albedo] = multisensor(90,(0.5:4.5)*360/5);
                elseif N == 2
                % sunto2: 21 sky regions (2 per sensor, except zenith), 10 albedo wedges
                    [~,obj.sky] = multisensor(0,0,45,(0:9)*360/10,90,(0:9)*360/10);
                    [~,obj.albedo] = multisensor(90,(0:9)*360/10);
                    obj.info.N = N;
                    obj.info.name = 'sunto2';
                else, error(errmsg);
                end
            case 'bucky'
            % 71 semi-regular sky regions, 6 dodecahedral albedo regions
                assert(isempty(N),errmsg);
                obj.info.name = 'bucky';
                a = atand(2);
                w = (0:4)*72;
                n12 = multisensor(0,0,a,w,180-a,w+36,180,0);    % starting with an icosaedron...
                n = [n12,voronoisphere(n12,'resolution',Inf)];  % add vertices of dual dodecahedron
                n = [n,voronoisphere(n,'resolution',Inf)];      % refine once more

                n(:,n(3,:)<sind(15)) = [];
                n = [n,multisensor(90,(0:19)*360/20,75,(0.5:19.5)*360/20)];
                [~,~,P] = voronoisphere(n,'resolution',1.5);
                obj.sky = cellfun(@polygon3d,P','unif',0);

                n12(:,n12(3,:) > sind(15)) = [];
                [~,~,P] = voronoisphere(n12,'resolution',1.5);
                obj.albedo = cellfun(@polygon3d,P','unif',0);
            case 'tin'
            % Recursive division of the faces of an icosahedral tessellation, to achieve 40,
            % 160,640,2560,... triangular faces on the sky hemisphere.
            % Wedges (two triangles wide) are used for albedo

                N = log(N/10)/log(4);
                if mod(N,1) ~= 0
                    N = max(1,round(N));
                    warning(['tinN works only for integers N = 10·4^z, for z > 0. Using ',...
                     'the closest match tin%d. Use ur%d if not happy.'],10*4^N,N);
                end
                obj.info.N = N;
                obj.info.name = sprintf('tin%d',N);

                n = spherepoints(2 + 10*4^N,'regular',true);

                % Keep triangles on the upper hemisphere as sky-regions
                T = convhull(n);
                T(any(reshape(n(T,3) < 0,[],3),2),:) = [];
                obj.sky = arrayfun(@(j) polygon3d(n(T(j,:)',:)'),1:size(T,1),'unif',0)';

                % Use a wedge-division with half the points on the horizon, for albedo
                h = abs(n(:,3)) < eps(1);
                [~,idx] = sort(atan2(n(h,2),n(h,1)));
                h = n(h,:);
                h = h(idx(1:2:end),:);
                [~,~,P] = voronoisphere(h,'resolution',1);
                obj.albedo = cellfun(@polygon3d,P','unif',0);

            case 'rb'
            % Recursive bisection distribution: 'RBN' = RB(N)
            % 1. Staring at horizon, add 'square' elements every 360°/N, move up to r = 360°/N
            % 2. While zenith is not reached, keep adding elements every d = 360°·sind(r)/n 
            %    making n(i+1) = n(i)/2, and r(i+1) = r(i) + d(i). 
            %    If n(i) is not divisible by 2, stop.
            % 3. Create a voronoi tesellation from the centers of those elements.

                assert(N > 4 && mod(N,4) == 0,'RBXX requires 4·z integer with z > 1');
                obj.info.N = N;
                obj.info.name = sprintf('rb%d',N);

                [~,~,P] = voronoisphere(multisensor(90,(0:N/2-1)*360/(N/2)),'resolution',1.5);
                obj.albedo = cellfun(@polygon3d,P','unif',0);

                r = 90; z = zeros(3,0); c = []; d = 0; n = N;
                for j = 1:ceil(log2(N/4))+1
                    d = max(d,360*sind(r)/n);
                    if round(r/d) > 1 && mod(n,2) == 0
                        z = [z,multisensor(r-d/3,(0:n-1)*360/n)]; %#ok
                        r = r - d;
                        n = n/2;
                        c = [c,round(d/2,1)]; %#ok
                    else
                        n = round(360*sind(r)/d);
                        if n > 1, z = [z,multisensor(2*r/3,(0:n-1)*360/n)]; %#ok
                        else, z = [z,[0;0;1]]; %#ok
                        end
                        break;
                    end
                end
                [~,~,P] = voronoisphere(z,'resolution',1.5);
                obj.sky = cellfun(@polygon3d,P','unif',0);
            case 'ur'
                obj.info.N = N;
                obj.info.name = sprintf('ur%d',N);
                X = spherepoints(N,'symmetric',true);
                [~,~,P] = voronoisphere(X,'resolution',1.5);
                obj.sky = cellfun(@polygon3d,P','unif',0);

                B = X(cellfun(@(p) any(p.z < 0),obj.sky),:);
                B(:,3) = 0;
                [~,~,P] = voronoisphere(B,'resolution',1.5);
                obj.albedo = cellfun(@polygon3d,P','unif',0);
            otherwise
                error(errmsg);
            end 
        end
        
        function x = rings(R,N,tol)
        % Return a cell-array of sum(N) 3D-polygons, resulting from dividing each spherical segment
        % R(j) < sqrt(x²+y²) < R(j+1) into N(j) segments. Arc segments have a maximum length TOL.
        %
        % EXAMPLE: for an orthographic projection of the 13 Satel-Light zones:
        %   C = RINGS([0 5/13 sqrt(105)/13 1],[1 4 8],angtol*pi/180)

            narginchk(3,3);
            assert(numel(R) == numel(N)+1,'Expecting numel(N)+1 vector for R');
            assert(all(diff(R)>0),'Expecting strictly sorted vector R');

            x = cell(numel(N),1);
            ns = zeros(numel(N),1);

            for j = 1:numel(N)
                if j == numel(N), ns(j) = N(j);
                else, ns(j) = lcm(N(j),N(j+1));
                end
                ns(j) = ceil(max(ns(j),2*pi*R(j+1)/tol)/ns(j))*ns(j);

                th = (0:1:ns(j)/N(j))/ns(j)*360;
                r = repmat(R(j+1),1,numel(th));
                if R(j) > 0
                % add return (inner) arc
                    th = [th,(ns(j-1)/N(j):-1:0)/ns(j-1)*360]; %#ok<AGROW>
                    r(end+1:numel(th)) = R(j);
                elseif N(j) > 1
                % add center vertex
                    th = [th,0];  %#ok<AGROW>
                    r(end+1) = 0; %#ok<AGROW>
                end
                z = sqrt(1-r.^2);
                x{j} = arrayfun(@(th0) polygon3d(th+th0,r,z,'pol'),(0:N(j)-1)/N(j)*360,'unif',0);
            end
            x = [x{:}];
        end
    end
    
    function obj = rotate(obj,z)
    % B = ROTATE(A,Z) - rotate static regions of A by z degrees around zenith.
    %
    % NOTE: this is not meant to be used for surface- or point-of-view dependent calculations,
    %   but to set the absolute orientation (relative to geographic North) of the sky regions.
    %   Use OBJ.PROJECTSTATIC and OBJ.PROJECTSOLAR for 3D POV rotations.
    
        validateattributes(z,{'numeric'},{'real','finite','scalar'},'','angle');
        R = polygon3d.rotmat(z,'Z');
        obj.sky = cellfun(@(p) polyrotate(p,R),obj.sky,'unif',0);
        obj.albedo = cellfun(@(p) polyrotate(p,R),obj.albedo,'unif',0);
        if obj.hasimplicit.sky
            obj.implicit.sky = cellfun(@(p) polyrotate(p,R),obj.implicit.sky,'unif',0);
        end
        if obj.hasimplicit.albedo
            obj.implicit.albedo = cellfun(@(p) polyrotate(p,R),obj.implicit.albedo,'unif',0);
        end
    end
    
    function [prj,rot,prj0] = projectstatic(obj,PolyPrj,horizon,R)
    % [PRJ,ROT,PRJ0] = PROJECTSTATIC(0BJ,POLYPRJ,HOR,R) - Project and clip sky/albedo regions
    %   in SHADINGREGIONS object OBJ, using the projection POLYPRJ, for a rotation of the
    %   view-pont R, and a local horizon profile HOR.
    %
    % INPUT: 
    %   OBJ - SHADINGREGIONS object
    %   POLYPRJ - POLYPROJECTOR object
    %   HOR - POLYGON3D local horizon profile
    %   R - 3·3 rotation matrix (OBJ region coordinates will be left-multiply by R).
    %
    % OUTPUT: 
    %   PRJ - Structure with fields {sky,albedo,horizon}. The first two are cell-arrays of
    %       POLYGON objects with the (flat) projected sky and albedo regions. The third is the
    %       flat POLYGON projection of HOR (polygon containing the visible sky).
    %   ROT - Structure with fields {sky,albedo,horizon}, contains the rotated and refined
    %       3DPOLYGON objects, which are an intermediate result of the projection.
    %   PRJ0 - Structure with fields {sky,albedo}, projected POLYGON regions, prior to clipping
    %       with projected horizon profile.
    %
    % See also: PROJECTSOLAR, POLYPROJECTOR, MOUNTROTATIONS

        if nargin < 3 || isempty(horizon) 
            horizon = polygon3d(0:PolyPrj.angtol/2:360,0,1,'sph'); % flat horizon, if missing
        end
        if nargin < 4 || isempty(R), R = eye(3); end
        
        if ~isequal(R,eye(3))
            obj.sky = cellfun(@(p) polyrotate(p,R'),obj.sky,'unif',0);
            obj.albedo = cellfun(@(p) polyrotate(p,R'),obj.albedo,'unif',0);
            horizon = polyrotate(horizon,R');
        end
        z = R(3,:);
        %flip = sign(R(3,3));
        
        % rotated viewlim is just viewlim
        viewlim = polygon(0:PolyPrj.angtol/2:360,PolyPrj.r0*PolyPrj.fun(PolyPrj.cX),'pol');
                
        % Generate and project a sample horizon-profile (polygon containing the sky)
        [prj.horizon,rot.horizon] = PolyPrj.project(horizon,z,true);
        prj_gnd = substractpolygons(viewlim,prj.horizon);
        
        batchproject = @(f) arrayfun(...
            @(j) PolyPrj.project(obj.(f){j},z,obj.haszenith.(f)(j)),1:numel(obj.(f)),'unif',0);
        
        % Project and clip sky/albedo regions
        [prj0.sky,rot.sky] = batchproject('sky');
        prj.sky = cellfun(@(p) intersectpolygons(p,prj.horizon),prj0.sky,'unif',0);
        if obj.hasimplicit.sky
            p = [prj0.sky{:}];
            if isempty(p)
                prj.sky{end+1} = prj.horizon;
                if nargout > 2, prj0.sky{end+1} = viewlim; end
            else
                prj.sky{end+1} = substractpolygons(prj.horizon,p,'pos','pos');
                if nargout > 2, prj0.sky{end+1} = substractpolygons(viewlim,p,'pos','pos'); end
            end
        end
        
        [prj0.albedo,rot.albedo] = batchproject('albedo');
        prj.albedo = cellfun(@(p) substractpolygons(p,prj.horizon),prj0.albedo,'unif',0);
        if obj.hasimplicit.albedo
            p = [prj0.albedo{:}]; 
            if isempty(p)
                prj.albedo{end+1} = prj_gnd;
                if nargout > 2, prj0.albedo{end+1} = viewlim; end
            else
                prj.albedo{end+1} = substractpolygons(prj_gnd,p,'pos','pos');
                if nargout > 2, prj0.albedo{end+1} = substractpolygons(viewlim,p,'pos','pos'); end
            end
        end
    end
    
    function [prj,rot] = projectsolar(obj,PolyPrj,s,R,horizon)
    % [PRJ,ROT] = PROJECTSOLAR(0BJ,POLYPRJ,S,R,HOR) - Project and clip circumsolar regions in
    %   SHADINGREGIONS object OBJ, using the projection POLYPRJ, for solar position vector S,
    %   a rotation of the view-pont R, and a local horizon profile HOR.
    %
    % INPUT: 
    %   OBJ - SHADINGREGIONS object
    %   POLYPRJ - POLYPROJECTOR object
    %   S - 3·N - unit vector(s) pointing towards the sun
    %   R - 3·3 - view-point rotation matrix (OBJ region coordinates will be left-multiply by R),
    %       or 3·3·M set of POV rotation matrices.
    %   HOR - POLYGON3D local horizon profile
    %
    % OUTPUT: 
    %   PRJ - {C,N,M} cell-array of POLYGON objects with the (flat) projected circumsolar regions,
    %       where C = obj.n.solar, and N, M is the number of sun-vectors & rotation matrices.
    %   ROT - {C,N,M} cell-array of rotated and refined 3DPOLYGON objects. Rotated objects are
    %       usually an intermediate result of the projection, but can represent overhead when the
    %       projection is outside the clipping area of PRJ.
    %
    % See also: PROJECTSTATIC, POLYPROJECTOR, MOUNTROTATIONS

        narginchk(2,5);
        if isempty(obj.cs), prj = {}; rot = {}; return; end
        
        if nargin < 3 || isempty(s), s = [0;0;1]; end
        if nargin < 4 || isempty(R), R = eye(3); end
        assert(size(s,1) == 3,'S must be a [3 Ns Nr] array');
        
        Nr = size(R,3);
        Ns = size(s,2);
        Nc = obj.n.solar;
        
        % PROVISIONAL: will this ever happen??
        if obj.hasimplicit.solar
            obj.solar = [{obj.implicit.solar},obj.solar];
            obj.cs = [16/60,obj.cs];
        end

        % visible(i,j,k) - CS disc i will not show at sun-position j for rotation k
        cosIA = permute(R(:,3,:),[3 1 2])*s; % [Nr Ns]
        visible = acosd(permute(cosIA,[3 2 1]))-obj.cs(:) < acosd(PolyPrj.cX); % [Nc Ns Nr]
        
        if ~any(visible) && nargout < 2
           prj = repmat({polygon.empty},Nc,Ns,Nr); return;
        end
        
        % Maximum horizon height
        if nargin < 5 || isempty(horizon); maxhz = 0; horizon = [];
        else, maxhz = max(horizon.z);
        end
        
        % Minimum z of visible (i.e. projected) point, for each rotation R
        minprjz = cos(acos(R(3,3,:)) + acos(PolyPrj.cX)); % [1 1 Nr]
        
        % clipped(i,j,k) - CS disc i might be clipped at sun-position j for rotation k
        clipped = visible & s(3,:) < sind(asind(maxhz) + obj.cs(:)); % [Nc Ns Nr]
        clipped = clipped & minprjz < maxhz;  % horizon is visible, otherwise prj does the work
       
        if isempty(horizon) && any(clipped(:))
            horizon = polygon3d(0:PolyPrj.angtol/2:360,0,1,'sph'); % flat horizon
        end

        prj = cell(Nc,Ns,Nr);
        rot = cell(Nc,Ns,Nr);
        [prj{~visible}] = deal(polygon.empty);

        for k = 1:Nr
            if nargout < 2 && ~any(visible(:,:,k),'all'), continue; end
            
            if any(clipped(:,:,k),'all')
                h = polyrotate(horizon,R(:,:,k)');
                prj_horizon = PolyPrj.project(h,R(3,:,k)',true);
            end 

            for j = 1:Ns
                if nargout < 2 && ~any(visible(:,j,k)), continue; end
                
                Rs = R(:,:,k)'*[null(s(:,j)'),s(:,j)]; % rotation for solar discs

                for i = Nc:-1:1 % start from largest (*)
                    rot{i,j,k} = polyrotate(obj.solar{i},Rs);
                    if ~visible(i,j,k), continue; end
                    
                    [prj{i,j,k},rot{i,j,k}] = PolyPrj.project(rot{i,j,k},Rs(:,3),true);
                    if ~clipped(i,j,k), continue; end
                    
                    if i == Nc
                    % clip with projected horizon
                        prj{i,j,k} = intersectpolygons(prj{i,j,k},prj_horizon);
                    else
                    % (*) use (hopefully simpler) bigger CS disc as mask
                        prj{i,j,k} = intersectpolygons(prj{i,j,k},prj{i+1,j,k});
                    end
                    if isempty(prj{i,j,k})
                    % (*) if CS disc i is not visible, neither are 1..i-1
                        if nargout > 1
                        % let loop continue to get missing rot{:,j,k}
                            visible(1:i-1,j,k) = false;
                        else
                            [prj{1:i-1,j,k}] = deal(polygon.empty);
                            break;
                        end
                    end
                end
            end
        end
    end
    
    function varargout = plot(obj,varargin)
    % OBJ.PLOT() - Lambert-Azimuthal-Equal-Area projection plot of SHADINGREGIONS object.
    % OPT = OBJ.PLOT() - Get the structure of default plotting options after plotting.
    % OPT = OBJ.PLOT(OPT) - Pass a (partial) structure of plotting options.
    % OPT = OBJ.PLOT(..,NAME,VAL..) - (partial) list of plotting options as name-value pairs.
    %
    %   Recognized options [default values] are:
    %
    %     OPT.rotation [eye(3)] - (starting) view point rotation matrix
    %     OPT.dynamic [false] - allow UI rotation with arrows after plotting
    %     OPT.horizon [random example horizon] 
    %     OPT.sunpos [0,0] - [az,el] (degrees), or [3,1] sun-position vector 
    %     OPT.colors.sky [blueish randomized colors], [OBJ.n.sky,3] colors for sky regions
    %     OPT.colors.gnd [greenish .. colors], [OBJ.n.albedo,3] colors for ground regions
    %     OPT.colors.sun [yellowish .. colors, 50% alpha], [OBJ.n.solar,4] colors for CS regions
    %     OPT.colors.edges ['none'] polygon edge colors
    %     OPT.figh [GUIfigure('ShadingRegions')]
    %     OPT.labels [false] plot numeric labels for each visible region (slow)
    %     OPT.wires [false] plot 3D wire mesh
    %     OPT.prj [Lambert projection], POLYPROJECTOR object.
    %
    %   Boolean propperties {dynamic, wires, labels}, work also as -flags, e.g. '-dynamic'.
        
        if isempty(varargin), varargin = {struct()}; end
 
        OPT.rotation = eye(3);
        OPT.dynamic = false;
        OPT.horizon = []; % samplehorizon();
        OPT.sunpos = [0,0];
        OPT.colors.sky = randomcolormap([0.5 0.7 1],obj.n.sky,0.2);
        OPT.colors.gnd = randomcolormap([0.6 0.8 0],obj.n.albedo,0.2);
        OPT.colors.sun = randomcolormap([1 1 0],obj.n.solar,0.2);
        OPT.colors.sun(:,4) = 0.5;
        OPT.colors.edges = 'w';
        OPT.figh = [];
        OPT.labels = false;
        OPT.wires = false;
        OPT.prj = 'lambert';        
        %OPT.prj = polyprojector('IAM');
        
        if numel(varargin) == 1 && isstruct(varargin{1})
            OPT = completestruct(varargin{1},OPT);
        elseif ~isempty(varargin)
            [bopt,varargin] = getflagoptions(varargin,{'-dynamic','-wires','-labels'});
            OPT = getpairedoptions(varargin,completeoptions(bopt,OPT),'restchk');
        end
        
        if numel(OPT.sunpos) == 2
            OPT.sunpos = sph2cartV(OPT.sunpos(1),OPT.sunpos(2))';
        end
        
        if ~isa(OPT.prj,'polyprojector'), OPT.prj = polyprojector(OPT.prj); end
        
        if isempty(OPT.figh) || ~ishandle(OPT.figh), OPT.figh = GUIfigure('ShadingRegions'); end

        % Project and clip sky/albedo regions
        [prj,rot] = projectstatic(obj,OPT.prj,OPT.horizon,OPT.rotation);

        % Rotate, project and clip circumsolar regions
        [prj.solar,rot.solar] = projectsolar(obj,OPT.prj,OPT.sunpos,OPT.rotation,OPT.horizon);
        
        % Remove splinter polygons
        prj.sky(cellfun(@area,prj.sky)/(pi*OPT.prj.r0^2) < 1e-4) = {polygon.empty};
        prj.albedo(cellfun(@area,prj.albedo)/(pi*OPT.prj.r0^2) < 1e-4) = {polygon.empty};
        prj.solar(cellfun(@area,prj.solar)/(pi*OPT.prj.r0^2) < 1e-4) = {polygon.empty};

        axes(OPT.figh); cla(); hold on;
        h.sky = batchplot(prj.sky,OPT.colors.sky,OPT.colors.edges);
        h.albedo = batchplot(prj.albedo,OPT.colors.gnd,OPT.colors.edges);
        h.solar = batchplot(prj.solar,OPT.colors.sun,OPT.colors.edges);
        s = OPT.prj.prj(OPT.rotation'*OPT.sunpos);
        plot(s(1),s(2),'y+');
        %axis([-1.5,1.5,-1.5,1.5]);
        axis equal
        set(gca,'visible','off');
        
        if OPT.wires
            hold on
            wireplot = @(P,C) arrayfun(@(j) polyplot(P{j},'none',C(j,:)),1:numel(P),'unif',0);
            wireplot(rot.sky,OPT.colors.sky);
            wireplot(rot.albedo,OPT.colors.gnd);
            wireplot(rot.solar,OPT.colors.sun);
            %p = arrayfun(@(x,y) polyrotate(polygon3d(polygon(72)),[x,y],'XY'),[0 90 0],[0 0 90],'unif',0);
            %wireplot(p,[1 1 1]*0.8);
        end
        
        if OPT.labels
            L = labelpolygons(prj.solar,[],'color',mean(OPT.colors.sun,1)*0.6,'fontsize',8);
            L = labelpolygons(prj.sky,L,'color',mean(OPT.colors.sky,1)*0.6,'fontsize',8);
                labelpolygons(prj.albedo,L,'color',mean(OPT.colors.gnd,1)*0.6,'fontsize',8);
        end
        
        if nargout > 0, varargout = {OPT,h}; end
        
        if OPT.dynamic, set(OPT.figh,'WindowKeyPressFcn', @KeyPress); end
        
        function h = batchplot(P,C,e)
        % Plot all polygons P with colors C, edge color E
        
            if isempty(P), h = matlab.graphics.primitive.Patch.empty; return; end
            haveholes = cellfun(@(p) any([p.hole]),P);
            if any(haveholes)
                h{1} = arrayfun(@(k) polyplot(P{k},C(k,:),e),find(haveholes));
                h{2} = batchplot(P(~haveholes),C(~haveholes,:),e);
                h = cat(1,h{:});
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
            h = patch('Faces',F,'Vertices',V,'FaceVertexCData',C(idx,:),'EdgeColor',e,...
                    'FaceColor','flat','FaceAlpha',alpha,'LineWidth',0.1);
        end

        function KeyPress(~,evd)        
            DELTAX = -5;
            DELTAY = 5;
            x = 0; y = 0;
            switch evd.Key
                case 'leftarrow', x = -DELTAX;
                case 'rightarrow', x = DELTAX;
                case 'uparrow', y = DELTAY;
                case 'downarrow', y = -DELTAY;
                case 'space'
                    OPT.rotation = eye(3);
                    plot(obj);
                    return;
                otherwise, return; 
            end
            R = polygon3d.rotmat([x,y],'YX');
            OPT.rotation = OPT.rotation*R;
            plot(obj,OPT);
        end

%         W = cell(3,1);
%         h = rotate3d(OPT.figh);
%         h.ActionPreCallback = @plotmesh;
%         h.ActionPostCallback = @updateprj;
%         h.Enable = 'on';
%         
%         function plotmesh(~,evd)
%             OPT.rotation = (polygon3d.rotmat(evd.Axes.View,'ZX')*[1 0 0;0 0 1;0 -1 0])';
%             hold on
%             wireplot = @(P,varargin) cellfun(@(p) polyplot(p,'none',varargin{:}),P,'unif',0);
%             %wireplot = @(P,varargin) cellfun(@(p) plot3(p.x([1:end,1]),p.y([1:end,1]),p.z([1:end,1]),varargin{:}),P);
%             W{1} = wireplot(rot.sky,OPT.colors.sky);
%             W{2} = wireplot(rot.albedo,OPT.colors.gnd);
%             W{3} = wireplot(rot.solar,OPT.colors.sun);
%             p = arrayfun(@(x,y) polyrotate(polygon3d(polygon(72)),[x,y],'XY'),[0 90 0],[0 0 90],'unif',0);
%             W{4} = wireplot(p,[1 1 1]*0.8);
%         end
%         function updateprj(~,evd)
%             cellfun(@(w) delete([w{:}]),W);
%             OPT.rotation = polygon3d.rotmat(evd.Axes.View,'ZX')*[1 0 0;0 0 1;0 -1 0]*OPT.rotation;
%             plot(obj);
%         end

        function c = randomcolormap(center,n,var)
        % Return n random color variations around center color
            c = rand(n,3)*var - var/2;
            c = center + c - mean(c,1);
            c = min(max(c,0),1);
        end
        function L = labelpolygons(P,L0,varargin)
        % Find maximum-inscribed-circle positions to place labels for a set of polygons
        
            if nargin < 2 || isempty(L0), L0 = zeros(0,2); end
            L = L0;
            for k = 1:numel(P)
                c = maxinscribedcircle(P{k},L);
                if ~isempty(c)
                    L = cat(1,L,c);
                    text(c(1),c(2),num2str(k),varargin{:});
                end
            end
            % arrayfun(@(k) text(L(k,1),L(k,2),num2str(k),varargin{:}),1:size(L,1));
            %L = [L;L0];
        end
    end
    
    function varargout = description(obj)
        
        ff = {'sky','albedo','solar'};
        tags = {'sky-region','albedo-region','circumsolar-ring'};
        msg = cell(3,1);
        for j = 1:3
            f = ff{j};
            switch obj.n.(f)
            case 0, msg{j} = sprintf('no %ss',tags{j}); 
            case 1
                if obj.hasimplicit.(f)
                    msg{j} = sprintf('1 implicit %s',tags{j});
                else
                    msg{j} = sprintf('1 %s',tags{j});
                end
            otherwise
                if obj.hasimplicit.(f)
                    msg{j} = sprintf('%d %ss (1 implicit)',obj.n.(f),tags{j});
                else
                    msg{j} = sprintf('%d %ss',obj.n.(f),tags{j});
                end
            end
        end
        msg(cellfun(@isempty,msg)) = [];
        msg = strjoin(msg,', ');
        if ~isempty(obj.info) && isfield(obj.info,'name')
            msg = [obj.info.name,': ',msg,'.'];
        else
            msg = [upper(msg(1)),msg(2:end),'.'];
        end      
        if nargout > 0, varargout{1} = msg; else, fprintf('%s\n',msg); end
    end
    
    function [sky,gnd,sun] = viewfactors(SR,prj,n,s,horizon,varargin)
    % [SKY,GND,SUN] = VIEWFACTORS(SR,PRJ,N,S,HOR) - Calculate view-factors for geometry SR,
    %   based on projection PRJ, for surface normal(s) N and solar-position(s) S, and considering
    %   a single (static) horizon profile HOR.
    %
    %   View factors are defined as the projected area of a given region (onto a given surface 
    %   plane, and for a given solar position), normalized by the region's solid angle. 
    %   E.g. the view-factor SKY(t,m,j) for the static sky-region SR.sky{j}, corresponding to 
    %   a surface with normal N(:,t,m) is given by:
    %
    %       SKY(t,m,j) ~ P(t,m).sky(j)/SR.solidangles.sky(j), where
    %       P(t,m) = PROJECTSTATIC(SR,PRJ,HOR,N(:,t,m))
    %
    %   VIEWFACTORS is designed to work with large sets of surface-normals and solar-positions,
    %   using SUNPOSMESH to generate an interpolation grid and reduce the number of required
    %   projections. Use VIEWFACTORS(..,'-exact') to override this behavior, and perform (slow!)
    %   step-by-step calculations for every surface-normal & solar-position combination.
    %
    % [...] = VIEWFACTORS(..,'opt',value) - Pass optional arguments to SUNPOSMESH.
    %   Default 'reduce' for SUNPOSMESH is set to TRUE.
    %
    % [...] = VIEWFACTORS(..,'-verbose') - print status / details. 
    %
    % INPUT:
    %   SR - SHADINGREGIONS object
    %   PRJ - projection key (e.g. 'IAM','ortho',..) or POLYPROJECTOR object
    %   N - [3 Nt* Nm] array of surface normals (pointing away from surface), a set of Nm static
    %       surfaces can be represented as a [3 1 Nm] array.
    %   S - [3 Nt] array of sun-vectors (pointing towards the sun), for Nt solar positions.
    %   HOR - POLYGON3D object (cartesian coordinates on unit sphere surface).
    %
    % OUTPUT:
    %   SKY - [Nt*,Nm,OBJ.n.sky] (single) array of view-factors for static sky-regions. 
    %       SKY(t*,m,j) is the projected area of region OBJ.sky{j} for surface normal N(:,t*,m), 
    %       divided by OBJ.solidangles.sky{j}
    %   GND - [Nt,Nm,OBJ.n.albedo] array of view factors for static albedo regions. GND(t,m,j) is 
    %       the projected area of region OBJ.albedo{j} for surface normal N(:,t,m), normalized by 
    %       the projected area of the complete visible hemisphere.
    %   SUN - [Nt,Nm,OBJ.n.solar] array of view factors. SUN(t,m,j) is the projected area of the
    %       circumsolar ring OBJ.cs{j} for surface normal N(:,t,m) and sun position S(:,t).
    %       Circumsolar view factors are normalized by the projected area of their respective
    %       regions for a surface pointing directly at the sun.
    %
    % See also: POLYPROJECTOR, SHADINGREGIONS, POAIRRADIANCE, PROJECTSTATIC
    
        [opt,varargin] = getflagoptions(varargin,{'-verbose','-exact'});
        if opt.exact && ~isempty(varargin)
            warning('Optional arguments for SUNPOSMESH will be ignored, -exact flag set'); 
        end
        [meshopt,varargin] = getpairedoptions(varargin,{'reduce'},{true});
        meshopt = [{'reduce',meshopt.reduce},varargin];
        
        printif = @(varargin) opt.verbose && fprintf(['\t' varargin{1} '\n'],varargin{2:end});
    
        if nargin < 2 || isempty(prj), prj = polyprojector('ortho'); end
        if ischar(prj)
           try prj = polyprojector(prj); end %#ok<TRYNC>
           assert(isa(prj,'polyprojector'),'PRJ: Expecting POLYPROJECTOR or valid character-key');
        end
        
        % Check the scaling of the projection, compare a projected 5° spherical cap vs its
        % theoretical solid angle. For ortho/lambert/IAM the usual is 0 <= 1-a/a0 < ~1e-2
        d = 1-area(prj.project(polygon3d(0:1:360,85,1,'sph')))/(2*pi*(1-sind(85)));
        if d < 0 || d > 1e-2
            warning('Near-zenith projection scale (%0.3f) seems off, check normalization',1-d);
        end
        
        if nargin < 3 || all(size(n)==0), n = [0;0;1]; end % n = zeros(3,1) is kept!
        if nargin < 4 || all(size(s)==0), s = [0;0;1]; end % id.
        if nargin < 5 || isempty(horizon), horizon = []; 
        else
            assert(isa(horizon,'polygon3d'),'HORIZON: Expecting POLYGON3D or empty');
        end

        validateattributes(s,{'numeric'},{'finite','2d','size',[3,NaN]},'','S');
        validateattributes(n,{'numeric'},{'finite','size',[3,NaN,NaN]},'','n');
        Nt = size(s,2);
        Nm = size(n,3);
        assert(any(size(n,2) == [1,Nt]),'size(N,2) is expected to be 1 or Nt');
        
        Nn = size(n,2)*Nm;     % No. of sun-position / surface-normal combinations
        % if size(n,2) == 1, n = repmat(n,1,Nt); end
        
        if ~opt.exact
        % Use an interpolation mesh that covers all surface normals
        % ('reduce' flag ensures original points are used if that's less work)
            printif('Getting Surface-Normal mesh...');
            [n0,Wn] = sunposmesh(reshape(n,3,Nn)',meshopt{:}); % [Gn,3] and [Nn,Gn] arrays
        else
           n0 = reshape(n,3,Nn)';
           Wn = speye(Nn);
        end
        Gn = size(n0,1);
        
        if Gn == Nn
            printif('... using original %d surface normals.',Nn);
        else
            printif('... %d mesh vertices for %d surface normals.',Gn,Nn);
        end
        
        % Static (sky/albedo) regions
        printif('Calculating static (%d-sky, %d-gnd) view-factors...',SR.n.sky,SR.n.albedo);
        F_sky = zeros(Gn,SR.n.sky);    
        F_gnd = zeros(Gn,SR.n.albedo);
        for k = 1:Gn
            % Project static regions at each mesh point
            R = [null(n0(k,:)),n0(k,:)'];
            p = projectstatic(SR,prj,horizon,R);
            F_sky(k,:) = cellfun(@area,p.sky)./SR.solidangles.sky;
            F_gnd(k,:) = cellfun(@area,p.albedo)./SR.solidangles.albedo;
        end
        F_sky = full(Wn*F_sky);   % ... and interpolate
        F_gnd = full(Wn*F_gnd);
        sky = reshape(single(F_sky),[],Nm,SR.n.sky);
        gnd = reshape(single(F_gnd),[],Nm,SR.n.albedo);
        
        sky = permute(sky,[1,3,2]);
        gnd = permute(gnd,[1,3,2]);
        
        % Circumsolar regions
        if ~isempty(SR.cs) && nargout > 2

            if ~opt.exact
            % Use an interpolation mesh that covers all solar-positions
                printif('Getting Sun-Position mesh...');
                [s0,Ws] = sunposmesh(s',meshopt{:}); % [Gs,3] and [Nt,Gs]
            else
                s0 = s';
                Ws = speye(Nt);
            end
            Gs = size(s0,1);
            
            if Gs == Nt
                printif('... using original %d sun-vectors.',Nt);
            else
                printif('... %d mesh vertices for %d sun-vectors.',Gs,Nt);
            end

            printif('Calculating circumsolar (%d) view-factors...',SR.n.solar);
            
            % Not all surface-normals are relevant for all sun-positions. Multiplying both
            % interpolation matrices produces non-zero values whenever a given surface-normal
            % grid-point is relevant for a given sun-position-grid-point:
            
            if Nn == Nm, Wn = repmat(Wn,Nt,1); end % static system
            coincident = full(spones(Wn)'*repmat(spones(Ws),Nm,1) > 0); % [Gn Gs]
            
            cosIA = n0*s0'; % [Gn Gs] array, cosine of incidence angle
            if isempty(horizon), maxhz = 0; else, maxhz = max(horizon.z); end

            % Horizon is visible: min. visible (i.e. projected) z < maxhz, for each normal n0
            vishz = cos(acos(n0(:,3)) + acos(prj.cX)) <  maxhz;
            
            % clipped(i,j,k) - CS disc i might be clipped at sun-position j for rotation k
            visible = acosd(cosIA)-shiftdim(SR.cs,-1) < acosd(prj.cX);     % [Gn Gs Ncs]
            clipped = s0(:,3)' < sind(asind(maxhz) + shiftdim(SR.cs,-1));  % [1 Gs Ncs]
            clipped = clipped & visible & vishz; 

            simple = ~any(clipped,3) & coincident;  % [Gn Gs]
            tricky = coincident & any(clipped,3);   % [Gn Gs]
            
            % A. Use a linear-interpolant for Incidence-Angle-driven combinations
            ia = linspace(0,90+SR.cs(end),ceil((SR.cs(end)+90)/prj.angtol))';
            if numel(ia) > nnz(simple)
                ia = acosd(cosIA(simple)); ia = ia(:);
                Wia = speye(numel(ia));
            else
                cias = cosIA(simple);
                Wia = interpmatrix(cias(:),cosd(ia));
            end
            
            printif('...Incident-Angle interpolation (%d prj)',numel(ia));
            p = projectsolar(SR,prj,[zeros(1,numel(ia));sind(ia)';cosd(ia)']); % {Nc,Nia}
            F_ia = cellfun(@area,p');
            F_ia(:,2:end) = diff(F_ia,1,2);
            F_ia = F_ia./SR.solidangles.solar;
            
            F_sun = zeros(Gn*Gs,SR.n.solar);
            F_sun(simple,:) = Wia*F_ia;
            F_sun = reshape(F_sun,Gn,Gs,SR.n.solar); % [Gn Gs Ncs]

            if any(tricky)
                % C. Do step-by-step projection for tricky cases
                printif('...Special-cases (%d prj)',nnz(tricky));
                for j = find(any(tricky,2)')
                    R = [null(n0(j,:)),n0(j,:)'];
                    p = projectsolar(SR,prj,s0(tricky(j,:),:)',R,horizon);
                    f = cellfun(@area,p');
                    f(:,2:end) = diff(f,1,2);
                    f = f./SR.solidangles.solar;
                    F_sun(j,tricky(j,:),:) = shiftdim(f,-1);
                end
            end
            
            % elev = repmat(atan2d(s0(:,3),rssq(s0(:,1:2),2)),1,Gn)';
            % scatter3(acosd(cosIA(:)),elev(:),reshape(F_sun(:,:,3),[],1),1,~tricky(:)); xlabel('IA'); ylabel('zenith');

            sun = zeros(Nt*Nm,SR.n.solar,'single');
            for j = 1:SR.n.solar
                x = Ws*F_sun(:,:,j)'; %[Nt Gn]
                for k = 1:Nm
                    block = (k-1)*Nt + (1:Nt)';
                    sun(block,j) = full(dot(Wn(block,:),x,2));
                end
            end
            sun = reshape(sun,Nt,Nm,SR.n.solar);
            sun = permute(sun,[1,3,2]);

            % cosIA = shiftdim(s(1,:,:).*n(1,:,:)+s(2,:,:).*n(2,:,:)+s(3,:,:).*n(3,:,:),1);
            % elev = repmat(atan2d(s(3,:),rssq(s(1:2,:),1)),Nm,1)';
            % tricky = reshape(dot(repmat(Ws*tricky',Nm,1),Wn,2),Nt,Nm);
            % scatter3(acosd(cosIA(:)),elev(:),reshape(sun(:,:,2),[],1),1,~tricky(:)); xlabel('IA'); ylabel('zenith');
        elseif nargout > 2
            sun = zeros(Nt,0,Nm,'single');
        end
    end
    
    function y = isequaltol(A,B,tol)
    % Y = ISEQUALTOL(A,B,TOL) - determine if objects A,B represent equal tessellations of the unit
    %   sphere within tolerance TOL (default SimOptions.RelTol).
    
        if nargin < 3, tol = getSimOption('RelTol'); end
        
        if ~isscalar(A) || ~isscalar(B)
            [A,B] = compatiblesize(A,B);
            y = arrayfun(@isequaltol,A,B,tol);
            return;
        end
        
        y = false;
        for f = {'hasimplicit','haszenith','n'}
            if ~isequal(A.(f{1}), B.(f{1})), return; end
        end
        
        isequalallright = @(x,y) all(abs((x-y)./(2*pi)) <= tol);
        if ~isempty(A.cs)
            if ~isequalallright(A.cs.^2,B.cs.^2), return; end
        end
        
        if ~isequalallright(A.solidangles.sky,B.solidangles.sky), return; end
        if ~isequalallright(A.solidangles.albedo,B.solidangles.albedo), return; end

        % Compare Lambert projection of sky/albedo regions
        prj = polyprojector('lambert');
        prjA = projectstatic(A,prj);
        prjB = projectstatic(B,prj);
        
        for f = {'sky','albedo'}
           a0 = max(cellfun(@area,prjA.(f{1})),cellfun(@area,prjB.(f{1})));
           ai = cellfun(@(a,b) area(intersectpolygons(a,b)),prjA.(f{1}),prjB.(f{1}));
           if ~all((a0-ai)/pi <= tol), return; end
        end

        y = true;
    end
end
methods (Static = true)
   function obj = loadobj(S)
   % OBJ = LOADOBJ(OBJ) - calculate and set all transient properties.
    
        isthere = @(s,f) isfield(s,f) | isprop(s,f);
        obj = ShadingRegions();
        if isthere(S,'info'), obj.info = S.info; end
        if isthere(S,'sky'), obj.sky = S.sky(:)'; end
        if isthere(S,'albedo'), obj.albedo = S.albedo(:)'; end
        if isthere(S,'cs'), obj.cs = sort(S.cs(:))'; end
        
        % Clip all sky-regions to horizon
        CPrj = polyprojector('azim','angtol',45,'clip',90);
        for j = 1:numel(obj.sky)
            [~,obj.sky{j}] = CPrj.project(obj.sky{j});
        end
    
        % Use a Lambert equal area projection to get approx. region solid angles
        AzPrj = polyprojector('lambert','normalize',sqrt(2));
        tol = (AzPrj.angtol*pi/180)^2/6; % area tolerance
        
        % Generate circumsolar-region concentric rings
        da = @(r) 360/ceil(360/min(0.8*AzPrj.angtol/sind(r),45));
        obj.solar = arrayfun(@(r) polygon3d(da(r)/2:da(r):360,90-r,1,'sph'),obj.cs,'unif',0);
        
        % Evaluate whether each region contains the projected zenith [0,0]
        for f = {'sky','solar','albedo'}
            Pj.(f{1}) = cellfun(@AzPrj.project,obj.(f{1}),'unif',0);
            aroundzenith = cellfun(@(p) insidepolygon(p,0,0),Pj.(f{1}));
            obj.haszenith.(f{1}) = aroundzenith & cellfun(@(p) ~any(p.z < 0),obj.(f{1}));
            if ~any(aroundzenith & ~obj.haszenith.(f{1})), continue; end
            for j = find(aroundzenith & ~obj.haszenith.(f{1}))
            % flip projection (make holes) for area checksum later on
                Pj.(f{1}){j} = AzPrj.project(obj.(f{1}){j},[0,0],false);
            end
        end
        assert(nnz(obj.haszenith.sky) <= 1,'At most 1-sky region should contain zenith');
        assert(~any(obj.haszenith.albedo),'No albedo regions should contain zenith');

        % Evaluate wheter the provided regions cover the complete dome within AzPrj.tol
        skydome = polygon(0:AzPrj.angtol/2:360,AzPrj.r0*AzPrj.fun(0),'pol');
        ground = polygon(0:AzPrj.angtol/2:360,AzPrj.r0*AzPrj.fun(-cosd(ShadingRegions.TINYARC)),'pol');
        % p = polygon3d(0:AzPrj.angtol/2:360,0,1,'sph');
        % skydome = AzPrj.project(p,[0,0],true);
        % ground = AzPrj.project(p,[0,0],false);
        
        obj.solidangles.sky = cellfun(@area,Pj.sky);
        obj.solidangles.albedo = cellfun(@area,Pj.albedo);
        
        obj.solidangles.solar = cellfun(@area,Pj.solar); % polygonal caps
        obj.solidangles.solar(2:end) = diff(obj.solidangles.solar);
        % A_sun = 2*pi*(1-cosd(obj.cs));               % area of spherical caps
        % A_sun(2:end) = A_sun(2:end)-A_sun(1:end-1);  % area of rings
        % obj.solidangles.solar = A_sun;

        % Flip scenario upside down for albedo regions
        % P = SR.projectstatic(prj,[],[1 0 0;0 -1 0;0 0 -1]);
        % A_gnd = cellfun(@(p) p.area,P.albedo);
        
        if ~isempty(obj.sky), skydome = substractpolygons(skydome,[Pj.sky{:}]); end
        if ~isempty(obj.albedo), ground = substractpolygons(ground,[Pj.albedo{:}]); end

        obj.hasimplicit.sky = area(skydome)/(2*pi) > tol;
        if obj.hasimplicit.sky
        % save the implicit missing region that fills the gaps
            obj.implicit.sky = AzPrj.inverse(skydome);
            obj.solidangles.sky(end+1) = area(skydome);
        else
            obj.implicit.sky = polygon3d.empty;
        end
        assert(abs(1-sum(obj.solidangles.sky)/(2*pi)) < tol,'Overlapping sky areas?');

        obj.hasimplicit.albedo = area(ground)/(4*pi) > 2*tol;
        if obj.hasimplicit.albedo
            obj.implicit.albedo = AzPrj.inverse(ground);
            obj.solidangles.albedo(end+1) = area(ground);
        else
            obj.implicit.albedo = polygon3d.empty;
        end
        assert(abs(1-sum(obj.solidangles.albedo)/(4*pi)) < 2*tol,'Overlapping ground areas?');
        
        % obj.hasimplicit.solar = true;
        % obj.implicit.solar = polygon3d(da(16/60)/2:da(16/60):360,90-16/60,1,'sph');
        obj.hasimplicit.solar = false;
        obj.implicit.solar = polygon3d.empty;
        
        obj.n.sky = numel(obj.sky) + obj.hasimplicit.sky;
        obj.n.albedo = numel(obj.albedo) + obj.hasimplicit.albedo;
        obj.n.solar = numel(obj.solar) + obj.hasimplicit.solar;
    end 
end
end











