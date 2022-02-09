classdef ShadingResults < matlab.mixin.Copyable
% SHADINGRESULTS - structure-class to hold the (packed) results of a Shading-Simulation.
%
% PROPERTIES:
%
% info - Structure with details about where the results come from, it should contain at least:
%   timing - STOPWATCH object with simulation-time info.
%   date - end-of-simulation
%   platform - result of VER('MATLAB')
%   code - result of GETGITINFO (structure with GIT branch and commit hash) 
%
% worldgeom - SHADINGREGIONS object, definining the geometry of the sky/albedo/circumsolar regions
% mountgeom - (packed) PVAREA object, describing mount geometry (modules, cell-blocks, cells)
%
% pts - copy of TRCK.analysedpoints: 3·Np array, diffuse shading evaluation points
% az, el - [Nt,1] vectors, solar positions for which analysis was performed.
%
% fullshaded - [Nt,Nu sparse bool] flags full-shading at time t, mount u
% partshaded - [Nt,Nu sparse bool] flags partial-shading.., Ns = nnz(partshaded)
%
% BLpoly - [Ni,1] List of beam-light polygons (i.e. complement of beam-shades over the mount plane)
% BLPidx - [Ns,Ni sparse bool] pointer-indices to beam-light polygons BLPoly. Ns = nnz(partshaded)
%   For a step-by-step simulation, BLPidx is a sparse boolean matrix [speye(Ni), if no two
%   polygons are equal], where BLPidx(s,i) = 1 means that the polygon for shading condition s is
%   to be found on i. After interpolation, [w,i] = find(BLPidx(s,:)) means that, for shading 
%   conditions s, polygons i should be assigned weights w.
%
% BshF -[Ni,Nm uint16] module linear-shading factors. If BLPidx(t,u) = k > 0, BshF(k,m)/(2^16 -1) 
%       yields the area fraction of module m of mount u covered by shade at t.
%
% Nsb - [Ni,Nm uint8] number of shaded cell-blocks on module m of mount u at t.
%
% DwF - Diffuse component weight factors, f(A) = projected-area(A) / solid-angle(A) · (2^16 - 1)
%       for each diffuse irradiance component bzw. OBJ.worldgeom region A, considering only the
%       effect of geometry (surface orientation), IAM, and far-horizon profile.
%       [Nu,1]** vector of structures, with fields {'sky','albedo','solar','gndbeam'}, each an
%       [Nt,N(c)]* UINT16 array, where N(c) = OBJ.worldgeom.n.(c) is the number of Diffuse   
%       Irradiance components (sky-regions) of each type. 
%
%       E.g. OBJ.DwF(u).sky(t,c) = 32768 implies that 32768/65535 ~ 25% of the radiance from 
%       sky-region c is incident on surface u at time-step t (in the absence of obstacles).
%
%       (*) For static mounts, DwF.sky and DwF.albedo are [1,N(c)] vectors.
%       (**) DwF is scalar (i.e. equal for all mounts) when all share a common surface normal
%       
% DshF - Shaded Diffuse component weight Factors. Equivalent to DwF, but considering the presence
%       of near obstacles (therefore Point-of-View-dependent)
%       [Nu,Npt] array of structures with fields {'sky','albedo','solar','gndbeam'} ,where Npt 
%       is the number of points analysed over each mount plane. Each field is an [Nt*,N(c)] UINT16 
%       array, where N(c) = OBJ.worldgeom.n.(c) is the number of Diffuse Irradiance components 
%       (sky-regions) of each type.
%
%       E.g. OBJ.DshF(u,p).sky(t,c) = 16384 implies that 16384/65535 ~ 25% of the radiance from 
%       sky-region c is incident on point p of surface u at time-step t, considering obstacles.
%
%       It is natural to expect that OBJ.DshF(u,p).(x)(t,c) <= OBJ.DwF(u).(x)(t,c) for all p,x
%       i.e. that near-shading decreases incidence irradiance. However, for albedo regions, the
%       presence of obstacles increases the reflective area, and the oposite might be the case.
%       Also, for OBJ.DwF(u).(x)(t,c)/(2^16 - 1) << 1, numerical precission issues can account 
%       for small apparent increases in DshF with respect to DwF. 
%
%       (*) For static mounts, DshF.sky and DshF.albedo are [1,N(c)] vectors. 
%
% DwF0 - Sensor diffuse component weight factors. [Ns,1] vector of structures, equivalent to
%       DwF but for Ns sensor orientations instead of Nu mount planes.
%
% The class is a pointer/handle to a structure, so beware when passing it around. Use B = COPY(A)
% and pass the copy, if you're not sure if the calling function will modify the variable.
%
% TODO: using arrays of structures OBJ.DshF(u,p).sky(t,c) makes for nicer notation, and in theory
%   was meant to allow for mounts with different numbers of analysis ponts (p). In praxis, it is
%   assumed that p is common, and (slow!) conversion to 4-D arrays is used to vectorize operations.
%   Should we just store everything as a 4D array with max(p)?
%
% TODO: clustering of mount "types" according to their direct shading "fingerprint" BshF(:,j)
%   should allow compression of the shading structure, and performance improvements in later
%   POA irradiance calculation.
%
% See also: SHADINGANALYSIS, POAIRRADIANCE, PVARRAYSOLVER, GUISHADING, GUIIRRTRANS, GUISOLVER
	
properties  
    info            % Structure with date/platform/version/timing
    worldgeom       % SHADINGREGIONS object
    mountgeom       % Packed PVAREA object
    pts             % [2,Np] coordinates of diffuse-shading analysis points on mount
    az              % [Nt·1 single] solar azimuth (degrees) 
    el              % [Nt·1 single] solar elevation (degrees)
    fullshaded      % [Nt,Nu sparse bool] flags full-shading at time t, mount u
    partshaded      % [Nt,Nu sparse bool] flags partial-shading.., Ns = nnz(partshaded)
    BLPidx          % [Ns,Ni sparse bool] pointer-indices to BLPoly
    BLPoly          % [Ni,1] Cell-array of 16-bit-packed polygons for partially-shaded mounts
    BshF            % [Ns,Nm] uint16 module linear-shading factors, Ns = nnz(partshaded)
    Nsb             % [Ns,Nm] uint8 number of shaded cell-blocks for each module
    DwF             % [Nu,1] structure, [Nt*·N(c) uint16] transposition-factor fields
    DwF0            % [Nk,1] structure, [Nt*·N(c) uint16] transposition-factors for sensors
    DshF            % [Nu,Np] structure, [Nt*×N(c) uint16] Diffuse shading factor fields
end
properties (Dependent = true)
    partidx         % [Nt,Nu sparse] returns row indices for BshF, Nsb, and BLPidx
    Nt              % scalar, number of simulation steps
    Nu              % scalar, number of mounts (trackers/tables)
    % Ns              % scalar, number of shaded conditions: nnz(partshaded)
    % Ni              % scalar, number of shaded polygons
    Np              % scalar, number of diffuse-shading evaluation points
    Nm              % scalar, number of modules per mount
    staticsystem
    equalrotations
end
methods
    function obj = ShadingResults(S)
    % OBJ = SHADINGRESULTS() - empty object constructor
    % OBJ = SHADINGRESULTS(S) - parse from existing structure or SHADINGRESULTS object S.
    %   Class constructor, meant for conversion of old result structures, or unpacked objects,
    %   e.g. the result of STRUCT(OBJ,'-unpack').
    %   For new results SHADINGANALYSIS should fill-up the structure with the propper format.
    %
    % See also: SHADINGANALYSIS, GUISHADING, SHADINGRESULTS.STRUCT
    
        UINT16MAX = 2^16 - 1;

        if nargin < 1 || isempty(S)
        % ShadingResults() - empty object constructor
            S.info = struct('timing',stopwatch(),...
                              'platform',system_dependent('getos'),...
                              'code',getGitInfo());
            %obj.worldgeom = ShadingRegions();
            S.az = zeros([0,1],'single');
            S.el = zeros([0,1],'single');
            S.DwF = struct('sky',{},'albedo',{},'solar',{});
            S.DwF0 = struct('sky',{},'albedo',{},'solar',{},'gndbeam',{});
            S.DshF = struct('sky',{},'albedo',{},'solar',{},'gndbeam',{});
            S.BLPidx = sparse(0,0);
            S.BLPoly = {};
            S.BshF = zeros([0,0],'uint16');
            S.Nsb = zeros([0,0],'uint8');
            return;
        end
        assert(isa(S,'ShadingResults') || isstruct(S),'Expecting structure or object');
        % ShadingResults(obj) - check & complete existing object
        
        isthere = @(s,f) isfield(s,f) | isprop(s,f);
        
        parsestruct(S,{'az','el'},'-n','-f','-r','-v','-e');
        parsestruct(S,{'BshF','Nsb'},'-n','-f','-r','-e');
        parsestruct(S,{'fullshaded','partshaded'},'-l','-e');
        
        % Enforce types
        if ~isa(S.az,'single'), S.az = single(S.az); end
        if ~isa(S.el,'single'), S.el = single(S.el); end
        if ~isa(S.BshF,'uint16')
            assert(all(S.BshF(:) >= 0 & S.BshF(:) <= 1),'BshF out of range');
            S.BshF = uint16(S.BshF*UINT16MAX); 
        end
        if ~isa(S.Nsb,'uint8'), S.Nsb = uint8(S.Nsb); end
        
        [Nt,Nu] = size(S.fullshaded);
        assert(numel(S.az) == Nt,'Non-matching az/el and fullshaded/partshaded');

%         % % (Re) calculate thresholds
%         % THRESH = getSimOption('cellshadingthreshold');
%         % THRESH = uint16(UINT16MAX*THRESH/prod(S.mountgeom.dims(2:end)));
%         THRESH = 0;
% 
%         if ~isthere(S,'fullshaded'), S.fullshaded = all(S.BshF >= UINT16MAX - THRESH,3); end
%         if ~isthere(S,'partshaded'), S.partshaded = any(S.BshF <= THRESH,3); end
%         assert(isequal(size(S.fullshaded),[Nt,Nu])        
        
        if ~isthere(S,'BLPidx') || isempty(S.BLPidx)
        % Build beam-shading-polygon indices from (unpacked) BLPoly

            if isequal(size(S.BLPoly),[Nt,Nu])
                S.BLPoly = S.BLPoly(S.partshaded);
            elseif nnz(S.partshaded) == 0
                S.BLPidx = speye(0);
                S.BLPoly = {};
            else
                assert(nnz(S.partshaded) == numel(S.BLPoly),'Non-matching BLPoly: missing index');
            end
            
            % [S.BLPoly,~,idx] = uniquecell(S.BLPoly);
            % S.BLPidx = sparse(1:numel(S.BLPoly),idx,1);
            S.BLPidx = speye(numel(S.BLPoly));
        else
            if ~issparse(S.BLPidx), S.BLPidx = sparse(S.BLPidx); end
            assert(size(S.BLPidx,1) == nnz(S.partshaded) && size(S.BLPidx,2) == numel(S.BLPoly),...
                'BLPidx must be an [Ns,Ni] sparse matrix');  
        end
        if ~issparse(S.partshaded), S.partshaded = sparse(S.partshaded); end
        
        % Reduce unpacked beam-shading-dependent propperties
        if size(S.BshF,1) ~= nnz(S.partshaded) || ~ismatrix(S.BshF)
            assert(size(S.BshF,1)*size(S.BshF,2) == Nt*Nu,'Non-matching BshF/partshaded');
            S.BshF = reshape(S.BshF,Nt*Nu,[]);
            S.BshF = S.BshF(S.partshaded(:),:);
        end
        if size(S.Nsb,1) ~= nnz(S.partshaded) || ~ismatrix(S.Nsb)
            assert(size(S.Nsb,1)*size(S.Nsb,2) == Nt*Nu,'Non-matching Nsb/partshaded');
            S.Nsb = reshape(S.Nsb,Nt*Nu,[]);
            S.Nsb = S.Nsb(S.partshaded(:),:);
        end

        usedpoly = any(S.BLPidx ~= 0,1);
        if ~all(usedpoly)
        % Remove unused polygons
            S.BLPoly = S.BLPoly(usedpoly);
            S.BLPidx = S.BLPidx(:,usedpoly);
        end
        
        % Parse/pack structures
        S.DwF = checkdiffusefactor(S.DwF,Nt,Nu);
        S.DwF0 = checkdiffusefactor(S.DwF0,Nt);
        S.DshF = checkdiffusefactor(S.DshF,Nt,Nu);
            
        % Check that polygons are packed
        if isempty(S.BLPoly) || isa([S.BLPoly{1}.x],'int16')
        % Nothing to do
        elseif isa(S.BLPoly{1},'polygon')
        % Pack polygons
            S.BLPoly = cellfun(@pack16,S.BLPoly,'unif',0);
        elseif isa([S.BLPoly{1}.x],'int64')
        % Re-pack polygons
            for j = 1:numel(S.BLPoly)
                for k = 1:numel(S.BLPoly{j})
                    f = S.BLPoly{j}(k).scale/polygon.SCALE16;
                    S.BLPoly{j}(k).x = int16(S.BLPoly{j}(k).x/f);
                    S.BLPoly{j}(k).y = int16(S.BLPoly{j}(k).y/f); 
                end
            end
        end

        % Complete info structure
        info = struct('platform','','code',struct('branch','','hash',''));
        if isthere(S,'timing'), info.timing = S.timing;
        else, info.timing = stopwatch();
        end
        if isthere(S,'info')
            if ~isstruct(S.info), S.info = struct('other',S.info); end
            S.info = completestruct(S.info,info);
        else
            S.info = info;
        end
        
        % Copy everything relevant from Structure/ShadingResults object
        % obj = ShadingResults();
        % p = getfield(metaclass(obj),'PropertyList');
        % p([p.NonCopyable]|[p.Abstract]|[p.Constant]|[p.Dependent]) = [];
        % p = {p.Name};
        p = {'info','worldgeom','mountgeom','az','el','fullshaded','partshaded','BLPidx',...
            'BLPoly','BshF','Nsb','DwF','DwF0','DshF'};
        for j = 1:numel(p)
            if ~isthere(S,p{j}), continue; end
            obj.(p{j}) = S.(p{j});
        end
        
        function x = checkdiffusefactor(x,Nt,M)
        % Verify that x is a structure with fields {sky,albedo,solar,[gndbeam]}, and each is a
        % [Nt,Nc] uint16-packed [0-1] array.
        
            if isempty(x), return; end
            
            if nargin < 3, M = size(x,1); end
            assert(isstruct(x) && any(size(x,1) == [0,1,M]),...
                '%s is not a propper-sized structure',inputname(1));
            
            FLD = {'sky','albedo','solar','gndbeam'};
            assert(all(isfield(x,FLD(1:3))),'Missing required field');
            if ~isfield(x,'gndbeam') || all(cellfun(@isempty,{x.gndbeam})), FLD(end) = []; end
                        
            for l = 1:numel(FLD)
                fld = FLD{l};
                Nc = size(x(1).(fld),2);
                for i = 1:numel(x)
                    v = x(i).(fld);
                    assert(size(v,2) == Nc && any(size(v,1) == [1,Nt]),...
                                        '%s(%d).%s has wrong size',inputname(1),i,fld);
                    if ~isa(v,'uint16')
                        assert(all(v(:) >= 0 & v(:) <= 1),...
                                        '%s(%d).%s out of range',inputname(1),i,fld);
                        x(i).(fld) = uint16(v*UINT16MAX);
                    end
                end
            end
        end
    end

    % Dependent propperty access methods
    function n = get.Nt(obj), n = size(obj.fullshaded,1); end
    function n = get.Nu(obj), n = size(obj.fullshaded,2); end
    function n = get.Nm(obj), n = size(obj.BshF,2); end
    function n = get.Np(obj), n = size(obj.DshF,2); end
    % function n = get.Ns(obj), n = nnz(obj.partshaded); end
    % function n = get.Ni(obj), n = numel(obj.BLPoly); end
    function a = get.staticsystem(obj), a = size(obj.DwF(1).sky,1) == 1; end
    function a = get.equalrotations(obj), a = numel(obj.DwF) == 1; end
    
    function yn = isequal(A,B)
        yn = false;
        if ~isa(B,'ShadingResults'), return; end
        % PROP = metaclass(A).PropertyList;
        % PROP = setdiff({props(~[props.Dependent]).Name},{'info','BLPoly','BLPidx'});
        PROP = {'BshF','DshF','DwF','DwF0','Nsb','az','el','fullshaded','mountgeom','partshaded','worldgeom'};
        eq = cellfun(@(p) isequal(A.(p),B.(p)),PROP);
        if ~all(eq), return; end

        if ~isequal(A.BLPoly,B.BLPoly)
        % One polygon index might just be a permutation of the other
            [~,idxA] = sortrows(A.BLPidx');
            [~,idxB] = sortrows(B.BLPidx');

            yn = isequal(A.BLPoly(idxA),B.BLPoly(idxB)) && ...
                 isequal(A.BLPidx(:,idxA),B.BLPidx(:,idxB));
            if ~yn, return; end
        end

         yn = all(eq);
    end
    
    function idx = get.partidx(obj)
        [r,c] = find(obj.partshaded);
        idx = sparse(r,c,1:numel(r),obj.Nt,obj.Nu);
    end

    function F = mountfilter(S,trckidx)
    % F = MOUNTFILTER(S,trckidx) - Reduce SHADINGRESULTS object S to a subset of mounts F.
    %   Resulting object is a copy of the original, with F.Nu = numel(idx)/nnz(idx), depending
    %   on whether trckidx is a numerical or logical filter.

        F = copy(S);

        F.fullshaded = F.fullshaded(:,trckidx);
        F.partshaded = F.partshaded(:,trckidx);
        
        sidx = nonzeros(S.partidx(:,trckidx));      % Reduce shading-dependent
        F.BshF = F.BshF(sidx,:);
        F.Nsb = F.Nsb(sidx,:);
        F.BLPidx = F.BLPidx(sidx,:);
        
        usedBLP = any(F.BLPidx,1);                  % Reduce polygon index
        F.BLPoly = F.BLPoly(usedBLP);
        F.BLPidx = F.BLPidx(:,usedBLP);

        if ~isempty(F.DwF), F.DwF = F.DwF(trckidx); end        % [Nt×Nu×4 double]
        if ~isempty(F.DshF), F.DshF = F.DshF(trckidx,:); end   % [Nt×Nu×Np×4 double]
    end
    
    function S = struct(obj,varargin)
    % S = STRUCT(OBJ) - convert SHADINGRESULTS object OBJ to a structure.
    % S = STRUCT(..,'-unpack') - expand BshF and Nsb to form [Nt,Nu,Nm] arrays, and unpack
    %       UINT16 fields DwF, DwF0, and DshF.
    % S = STRUCT(..,'except',FF) - Skip field(s) FF
        
        % MATLAB seems too stupid to tell appart struct(R) and struct(..,'field',R)
        if ~isa(obj,'ShadingResults'), S = builtin('struct',obj,varargin{:}); return; end
    
        [opt,varargin] = getflagoptions(varargin,{'-unpack'});
        opt.except = {};
        opt.only = {};
        opt = getpairedoptions(varargin,opt,'restchk');
        
        if ~isempty(opt.only)
            p = opt.only; 
        else 
            % p = getfield(metaclass(obj),'PropertyList');
            % p = {p.Name};
            p = {'info','worldgeom','mountgeom','az','el','fullshaded','partshaded','BLPidx',...
                'BLPoly','BshF','Nsb','DwF','DwF0','DshF','partidx','Nt','Nu','Np','Nm',...
                'staticsystem','equalrotations'};
        end
        if ~isempty(opt.except), p = setdiff(p,opt.except); end
        
        S(numel(obj)) = struct(); S = reshape(S,size(obj));
        for i = 1:numel(obj)
            % Copy all properties to structure
            for j = 1:numel(p), S(i).(p{j}) = obj(i).(p{j}); end

            if ~opt.unpack, return; end

            UINT16MAX = 2^16 - 1;

            % Unpack and expand BshF & Nsb
            if isfield(S,'BshF')
                S(i).BshF = zeros(obj(i).Nt*obj(i).Nu,obj(i).Nm,'single');
                S(i).BshF(obj(i).fullshaded,:) = 1;
                S(i).BshF(obj(i).partshaded,:) = single(obj(i).BshF)/UINT16MAX;
                S(i).BshF = reshape(S(i).BshF,obj(i).Nt,obj(i).Nu,obj(i).Nm);
            end

            if isfield(S,'Nsb')
                S(i).Nsb = zeros(obj(i).Nt*obj(i).Nu,obj(i).Nm,'uint8');
                S(i).Nsb(obj(i).fullshaded,:) = obj(i).mountgeom.dims(2);
                S(i).Nsb(obj(i).partshaded,:) = obj(i).Nsb;
                S(i).Nsb = reshape(S(i).Nsb,obj(i).Nt,obj(i).Nu,obj(i).Nm);
            end

            % Unpack structures
            if isfield(S,'DwF'), S(i).DwF = unpackstruct(S(i).DwF); end
            if isfield(S,'DwF0'), S(i).DwF0 = unpackstruct(S(i).DwF0); end
            if isfield(S,'DshF'), S(i).DshF = unpackstruct(S(i).DshF); end
        end

        function x = unpackstruct(x)
            for k = 1:numel(x)
                for f = fieldnames(x)
                    x(k).(f{1}) = single(x(k).(f{1}))/UINT16MAX;
                end
            end
        end
    end
    
    function B = interpolate(A,az,el,maxarc,minwt)
    % Interpolate Shading-Results
    % TODO: improve performance!
    
        if nargin < 4 || isempty(maxarc), maxarc = getSimOption('angtol'); end
        if nargin < 5 || isempty(minwt), minwt = getSimOption('RelTol'); end
        
        assert(numel(az) == numel(el),'Non-matching AZ,EL');
        az = az(:); el = el(:);
        
        tol = max(maxarc*minwt,15/3660)*pi/180; % min. tol. = 1 arc-second 
        if A.Nt == numel(el) && all(solarposition.arcdist(A.el,A.az,el,az,1) < tol)
            B = copy(A);
            B.el = el;
            B.az = az;
            return;
        end
                
        V = sph2cartV(A.az,A.el);
        s = sph2cartV(az,el); 
        [W,invalid] = interpmatrix(s,V,'-sph','maxarc',maxarc*pi/180,'tol',tol);
        W(abs(W) < minwt) = 0;
        invalid = invalid | ~any(W,2);
        W(~invalid,:) = W(~invalid,:)./sum(W(~invalid,:),2);
        invalid = invalid & el > 0;
        B = filterstructure(A,W);
        B.fullshaded(el <= 0,:) = true;
        B.az = az;
        B.el = el;

        if any(invalid)
           warning('%d/%d points outside interpolation tolerance of %0.2f°',nnz(invalid),B.Nt,maxarc);
           B.az(invalid) = NaN;
           B.el(invalid) = NaN;
        end
    end
    
    function B = filterstructure(A,filter,varargin)
    % Overloaded FILTERSTRUCTURE for SHADINGRESULTS objects
    % TODO: improve performance for sparse matrix filters (i.e. few timesteps from a dense full-
    %   year grid) - maybe prefiltering by any(filter,1). Also, the use of [Nu,..] struct arrays
    %   with [Nt,Nc] fields instead of a simple [Nt,Nu,Nc,..] probably causes a big overhead on
    %   the original FILTERSTRUCTURE
    
        assert(isscalar(A),'FILTERSTRUCTURE only works for scalar objects');
    
        SCALAR = {'info','mountgeom','worldgeom'};
        INDEXED = {'BLPidx','BLPoly'};
        DEPENDENT = {'partidx','Nt','Nu','Np','Nm','staticsystem','equalrotations'};
        
        matrixfilter = size(filter,1) > 1 && size(filter,2) > 1;
        if matrixfilter
            assert(size(filter,2) == A.Nt,'Expecting M by %d matrix filter',A.Nt);
            INDEXED = [INDEXED,{'fullshaded','partshaded'}]; % will be recalculated
        elseif islogical(filter)
            assert(numel(filter) == A.Nt,'Expecting %d logical vector index',A.Nt);
        else
            assert(all(filter <= A.Nt & mod(filter,1) == 0),'Index out of bounds');
        end
        
        % Copy to structure, unpack indexed fields and uint16'd values
        % remove confusing scalar fields and indexed fields, that need special treatment
        B = struct(A,'-unpack','except',[SCALAR,INDEXED,DEPENDENT]);
        if matrixfilter
        % Allow integer values to (provisonally) be averaged continuously
            B.Nsb = single(B.Nsb);
            B.s = sph2cartV(B.az,B.el);
        end
        
        B = filterstructure(B,filter,varargin{:});  % filter
        B.BLPoly = A.BLPoly;
        
        % Handle indexed fields
        if matrixfilter
            Ni = numel(A.BLPoly);
            Nf = size(filter,1);
            
            % Add dummy full-/no-shade polygons, for points that are partly-fullshaded, or 
            % partly-not-shaded, i.e. points j with 0 < FILTER(j,:)·fullshaded(:,u)) < 1 and/or
            % 0 < FILTER(j,:)·(~fullshaded(:,u) & ~partshaded(:,u)) < 1
            B.BLPoly(Ni+1) = {pack16(polygon())};
            B.BLPoly(Ni+2) = {A.mountgeom.border};
            
            % Create an expanded [Nf·Nu × Ni+2] B.BLPidx
            B.BLPidx = sparse(Nf*A.Nu,Ni+2);
            for u = 1:A.Nu
%                 [r,c] = find(obj.partshaded);
%                 idx = sparse(r,c,1:numel(r),obj.Nt,obj.Nu);
        
                iA = nonzeros(A.partidx(:,u));
                iB = find(A.partshaded(:,u));
                notshaded = ~(A.fullshaded(:,u) | A.partshaded(:,u));
                n = numel(iA) + nnz(A.fullshaded(:,u)) + nnz(notshaded); %(§)

                % Expand A.BLPidx(iA,:) into a full Nt × Ni+2 matrix that allows FILTER*X ...
                [i,j,v] = find(A.BLPidx(iA,:));
                X = sparse(iB(i),j,v,A.Nt,Ni+2,n);
                
                % (§) Include full-shaded and not-shaded points using dummy polygons Ni+1:2
                X(A.fullshaded(:,u),Ni+1) = 1; %#ok<SPRIX>
                X(notshaded,Ni+2) = 1; %#ok<SPRIX>
                
                % % The same, but without sparse indexing (slow according to MATLAB)
                % % ... doesn't seem to make much of a difference.
                % na = numel(iA);
                % nF = nnz(A.fullshaded(:,u));
                % n = na + nF + nnz(notshaded); %(§)
                % i = zeros(n,1); j = zeros(n,1); v = ones(n,1);
                % 
                % [ix,j(1:na),v(1:na)] = find(A.BLPidx(iA,:));
                % i(1:na) = iB(ix);
                % i(na+1:na+nF) = find(A.fullshaded(:,u));
                % i(na+nF+1:n) = find(notshaded);
                % j(na+1:na+nF) = Ni+1;
                % j(na+nF+1:n) = Ni+2;
                % % v(na+1:end) = 1;
                % X = sparse(i,j,v,A.Nt,Ni+2);

                % Place the filtered index FILTER*X in the corresponding rows of B.BLPidx
                B.BLPidx((u-1)*Nf+(1:Nf),:) = filter*X;
            end
                        
%             W = repmat(filter,1,A.Nu);  % expand to Nx by Nt·Nu
%             W = W(:,A.partshaded);      % reduce columns to Ns shade-cases, to match BLPidx
%             B.BLPidx = W*A.BLPidx;      % and get new weight matrix (Nx by Np)

            % % (Re) calculate thresholds
            THRESH = getSimOption('cellshadingthreshold');
            THRESH = THRESH/prod(A.mountgeom.dims(2:3));
            THRESH = max(THRESH,eps(1));

            B.fullshaded = all(B.BshF >= 1 - THRESH,3);
            B.partshaded = any(B.BshF >= THRESH,3) & ~B.fullshaded;
            B.Nsb = int8(ceil(B.Nsb - THRESH*A.mountgeom.dims(2)));
            
            B.az = atan2d(B.s(:,2),B.s(:,1));
            B.el = atan2d(B.s(:,3),hypot(B.s(:,1),B.s(:,2)));

            B.BLPidx = B.BLPidx(B.partshaded(:),:);
        else
        % Use filtered A.partidx
            idx = filterstructure(A.partidx,filter,varargin{:});
            B.BLPidx = A.BLPidx(nonzeros(idx),:);
        end
       
        for f = SCALAR, B.(f{1}) = A.(f{1}); end    % put back scalar fields

        B = ShadingResults(B); % ... and back to SHADINGRESULTS (remake index, etc.)
    end
    
    function X = mergestructures(S,idx,varargin)
    % Overloaded MERGESTRUCTURES for SHADINGRESULTS objects

        EQUALFIELDS = {'mountgeom','worldgeom','Nu','Np','Nm'};
        INDEXEDFIELDS = {'BLPidx','BLPoly'};
        
        if isempty(S), X = ShadingResults(); return; end
        
        for f = EQUALFIELDS
        % Check that invariant fields are consistent
            assert(all(arrayfun(@(s) isequal(s.(f{1}),S(1).(f{1})),S(2:end))),...
                'Non-matching field: %s',f{1});
        end
        
        X = struct(S,'-unpack','except',[INDEXEDFIELDS,EQUALFIELDS]);

        % Update partshaded indices to work with blkdiag(S.BLPidx)
        Ns = nnz(X(1).partshaded);
        for j = 2:numel(S)
           X(j).partidx(X(j).partshaded) = S(j).partidx(S(j).partshaded) + Ns;
           Ns = Ns + nnz(X(j).partshaded);
        end
        
        % Merge structures
        X = mergestructures(X,idx,varargin{:});

        % Filter/reorder BLPidx according to the merged partidx
        X.BLPidx = blkdiag(S.BLPidx);
        X.BLPidx = X.BLPidx(nonzeros(X.partidx),:);
        X.BLPoly = cat(1,S.BLPoly);

        % ... add fields back
        for f = EQUALFIELDS, X.(f{1}) = S.(f{1}); end
        
        X = ShadingResults(X); % ... and pack to go
    end
end % methods
methods (Static = true)
    function S = loadobj(S)
    % Resolve backwards compatibiliy issues
        
        isthere = @(s,f) isfield(s,f) | isprop(s,f);
        UINT16MAX = 2^16 - 1;
        packU16 = @(x) uint16(x*UINT16MAX);

        if ~isstruct(S.DwF) && ~isempty(S.DwF)
        % Update from (deprecated) single-array result structure w. Perez model
        % DwF = [Nt,Nu,Nc] where c = {Isotropic, Albedo, Circumsolar, HB}
        % DshF = [Nt,Nu,Np,Nc]
        % DshF0 = [Nt,Nc] = DwF0/DwF0(no-horizon) for single horizontal sensor

            parsestruct(S,{'DwF','DshF0_u16','DshF_u16','BiaF_u16'},'-n');

            DshF0 = single(S.DshF0_u16)/UINT16MAX;
            BiaF = single(S.BiaF_u16)/UINT16MAX;

            S.worldgeom = ShadingRegions('perez');

            % NOTE: in most Perez et al. formulations, HB component is zero for
            % horizontal planes, even though it projects as 1-cos(10°)
            S.DwF0.sky = mean(DshF0(:,[1,4]),1).*[1,1-cosd(10)];

            % Albedo was actually calculated as 1 - ISO/ISO(flat-horizon) for projected
            % isotropic areas ISO. So it already reflects projected visible ground.
            S.DwF0.albedo = mean(DshF0(:,2),1);

            % DwF_CS = Rb0_poa/Rb0_flat -> Rb0_hor = DshF_CS·Rb0_poa/DwF_CS
            Rb0 = nanmean(BiaF./S.DwF(:,:,3),2);
            S.DwF0.solar = DshF0(:,3).*Rb0;
            clear DshF0

            DshF = S.DshF_u16;
            S.DshF_u16 = [];
            DwF = single(S.DwF);
            S.DwF = struct('sky',{},'albedo',{},'solar',{},'gndbeam',{}); 
            S.DshF = struct('sky',{},'albedo',{},'solar',{},'gndbeam',{});

            for j = size(DwF,2):-1:1   
                S.DwF(j,1).sky = permute(DwF(:,j,[1,4]),[1,3,2]);                            
                S.DwF(j,1).albedo = permute(DwF(:,j,2),[1,3,2]);
                S.DwF(j,1).solar = permute(DwF(:,j,3),[1,3,2]);

                for k = size(S.DshF,3):-1:1
                    S.DshF(j,k).sky = (1-permute(DshF(:,j,k,[1,4]),[1,4,2,3])).*S.DwF(j,1).sky;
                    S.DshF(j,k).albedo = (1-permute(DshF(:,j,k,2),[1,4,2,3])).*S.DwF(j,1).albedo;
                    S.DshF(j,k).solar = (1-permute(DshF(:,j,k,3),[1,4,2,3])).*S.DwF(j,1).solar;
                end 
            end
            clear DwF DshF

            % Remove redundant factors (static regions for static modules)
            if size(unique([S.DwF.sky],'rows'),1) == 1
                for j = 1:numel(S.DwF), S.DwF(j).sky = S.DwF(j).sky(1,:); end
            end
            if size(unique([S.DwF.albedo],'rows'),1) == 1
                for j = 1:numel(S.DwF), S.DwF(j).albedo = S.DwF(j).albedo(1,:); end
            end
            if size(unique([S.DshF.sky],'rows'),1) == 1
                for j = 1:numel(S.DshF), S.DshF(j).sky = S.DshF(j).sky(1,:); end
            end
            if size(unique([S.DshF.albedo],'rows'),1) == 1
                for j = 1:numel(S.DshF), S.DshF(j).albedo = S.DshF(j).albedo(1,:); end
            end
        else
        % unpack int16 structures, if they're used
            if isthere(S,'DshF_u16'), S.DshF = S.DshF_u16; S.DshF_u16 = [];end
            if isthere(S,'DwF_u16'), S.DwF = S.DwF_u16; S.DwF_u16 = []; end
            if isthere(S,'DwF0_u16'), S.DwF0 = S.DwF0_u16; S.DwF0_u16 = []; end
        end

        if ~isempty(S.DwF) && ~isa(S.DwF(1).sky,'uint16')
        % Update objects with old DwF normalization (visible-dome instead of solid angles, 
        %   shading-fractions isntead of shaded-view-factors)
        % (commit b330109 and earlier ~June/2019)

            IAMprj = polyprojector('IAM','normalize',false);
            ViewLim_area = polygon(0:IAMprj.angtol/2:360,IAMprj.r0*IAMprj.fun(0),'pol');
            ViewLim_area = ViewLim_area.area;
            CS_areas = cellfun(@(p) area(IAMprj.project(p)),S.worldgeom.solar);

            for j = 1:numel(S.DwF0)
                S.DwF0(j).sky = packU16(S.DwF0(j).sky*ViewLim_area./S.worldgeom.solidangles.sky);
                S.DwF0(j).albedo = packU16(S.DwF0(j).albedo*ViewLim_area./S.worldgeom.solidangles.albedo);
                S.DwF0(j).solar = packU16(S.DwF0(j).solar.*CS_areas./S.worldgeom.solidangles.solar);
            end
            for j = 1:numel(S.DwF)
                S.DwF(j).sky = packU16(S.DwF(j).sky*ViewLim_area./S.worldgeom.solidangles.sky);
                S.DwF(j).albedo = packU16(S.DwF(j).albedo*ViewLim_area./S.worldgeom.solidangles.albedo);
                S.DwF(j).solar = packU16(S.DwF(j).solar.*CS_areas./S.worldgeom.solidangles.solar);
            end
            for j = 1:numel(S.DshF)
                for k = size(S.DshF,3):-1:1
                    S.DshF(j,k).sky = (1-S.DshF(j,k).sky).*S.DwF(j,1).sky;
                    S.DshF(j,k).albedo = (1-S.DshF(j,k).albedo).*S.DwF(j,1).albedo;
                    S.DshF(j,k).solar = (1-S.DshF(j,k).solar).*S.DwF(j,1).solar;
                end
            end
        end
        
        if ~isfield(S.DshF,'gndbeam'), S.DshF(end).gndbeam = []; end
        if ~isfield(S.DwF,'gndbeam'), S.DwF(end).gndbeam = []; end

        % Add SCALE field to packed polygons (not included in previous versions)
        if ~isempty(S.BLPoly) && isstruct(S.BLPoly{1}) && ~isfield(S.BLPoly{1},'scale')
           for j = 1:numel(S.BLPoly), [S.BLPoly{j}.scale] = deal(polygon.SCALE16); end
        end
        
        % Name changes
        if isfield(S,'geom'), S.worldgeom = S.geom; end
                
        if isempty(S.az) && evalin('base','exist(''SunPos'')')
        % Try to get missing solar position data from base
            S.az = single(evalin('base','SunPos.Az'));
            if numel(S.az) ~= size(S.BLPidx,1)
                warning('Unknown ShRes.az/el - non matching SunPos!');
                S.az = zeros([0,1],'single');
            else
                S.el = single(evalin('base','SunPos.El'));
            end
        end
        
        if evalin('base','exist(''Trackers'')')
            if isempty(S.mountgeom)
            % Create a (packed) copy of Trackers.geom, with tighter fits to pv-active-surfaces
                trckgeom = evalin('base','Trackers.geom');
                S.mountgeom = copy(trckgeom);
                for j = 1:trckgeom.dims(1)
                    S.mountgeom.elements(j).border = areaenvelope(S.mountgeom.elements(j));
                end
                S.mountgeom.border = areaenvelope(S.mountgeom);
                S.mountgeom.border = intersectpolygons(trckgeom.border,S.mountgeom.border);
                S.mountgeom = struct(S.mountgeom,'pack16'); % pack16 all border polygons
            end
            if isempty(S.pts)
                S.pts = evalin('base','Trackers.analysedpoints');
            end
        end

        if ~isempty(S.mountgeom) && (isempty(S.BshF) || isempty(S.Nsb))
        % (Re)calculate shading fractions on modules & cell-blocks
            shf = 1-pvArea.shadingfactors(S.mountgeom,S.BLPoly,'-pack16','depth',3);   
            S.BshF = uint16(mean(shf,3)*UINT16MAX);
            S.Nsb = uint8(sum(shf > 0,3));
        end
        
        if ~isthere(S,'fullshaded') || size(S.BLPidx,2) ~= numel(S.BLPoly)
        % Switch from sparse-double BLPidx (-1 as flag for full-shading) to sparse-boolean indices
            S.fullshaded = sparse(S.BLPidx < 0);
            S.partshaded = sparse(S.BLPidx > 0);
            idx = S.BLPidx(S.partshaded);
            if ~isequal(idx',1:nnz(S.partshaded))
            % re-order rows of BshF, Nsb
               S.BshF = S.BshF(idx,:);
               S.Nsb = S.Nsb(idx,:);
               S.BLPoly = S.BLPoly(idx);
            end
            S.BLPidx = speye(nnz(S.partshaded));
        end
        
        S = ShadingResults(S); 
    end
end
methods (Hidden = true)
    function varargout = plot(obj)
    % Interactive plot of diffuse and direct view- and shading-factors.

        ShR = struct(obj,'-unpack');
        
        % Store { DwF, DshF, 1-BshF } in F as [Nc*,Nt,Nu,Np*,Nm*] arrays
        F = {ShR.DwF,ShR.DshF};
        FLD = fieldnames(F{1});
        n = numel(FLD);
        for k = 1:2
        % Turn F(u,p).x(t,c) into F.x(c,t,u,p) 
            stack = @(varargin) reshape(cat(ndims(varargin{1})+1,varargin{:}),...
                [size(varargin{1}),size(F{k})]);
            F{k} = struct2cell(F{k}); % [n,p,q] cell of [Nt,Nc] matrices
            F{k} = arrayfun(@(j) permute(stack(F{k}{j,:}),[2,1,3,4]),1:n,'unif',0);
            F{k} = cell2struct(F{k}(:),FLD);
        end
        F{3} = 1 - permute(ShR.BshF,[4 1 2 5 3]); % [1,Nt,Nu,1,Nm]
        
        % General layout
        GUIfigure('ShadingResults','Shading Results','big');
        clf();
        ax = arrayfun(@(j) subplot(2,2,j),1:4);
        arrayfun(@(j) caxis(ax(j),[0,1]),1:4);
        
        set(ax(1),'Visible',false);
        colorbar(ax(1));

        % For Diffuse shading, plot DwF and DshF as colors on ShR.worldgeom plots
        [~,H{1}] = plot(ShR.worldgeom,'labels',true,'ax',ax(3));
        text(-1.5,1.4,'DwF','fontweight','bold');

        [~,H{2}] = plot(ShR.worldgeom,'labels',true,'ax',ax(4));
        text(-1.5,1.4,'DshF','fontweight','bold');

        H{3} = plotbeamshading(ax(2));
        
        % Draw drop-down lists that adjust reducing functions interactively:
        extendedlist = @(n) [{'min','mean','median','mode','std','max'},...
                                                arrayfun(@num2str,1:n,'unif',0)];
        plotcontrols({'popupmenu'},...
            {'Solar Position:','Mount:','Analysis Point:','Module:'},...
            arrayfun(extendedlist,[0,ShR.Nu,ShR.Np,ShR.Nm],'unif',0),...
            {'mean'}, @update,'position',ax(1).Position.*[1 1 0.8 1]);
        
        if nargout > 0, varargout{1} = H; end
        
        function H = plotbeamshading(ax)
        % TODO: plot beam-shading polygons for particular time-steps
            
            set(ax,'visible',false);
            
            NDIST = 2;
            THRESH = 0.05;

            % Attempt to build a regular TIN grid, to use TRISURF
            s = sph2cartV(ShR.az,ShR.el);
            T = convhull([double(s);[0,0,-1]]);
            T(any(T > size(s,1),2),:) = [];
            djk = @(j,k) rssq(s(T(:,j),:) - s(T(:,k),:),2);
            D = [djk(1,2),djk(2,3),djk(1,3)];
            f = any(D > NDIST*median(D),2);
            T(f,:) = [];
            regular = nnz(f)/numel(f) < THRESH;

            % subplot(ax)
            if ~regular
                % a = solarposition.fixtoplusminus180(ShR.az-90);
                % H = scatter(a,ShR.el,20,mean(ShR.BshF,2:3),'filled'); 
                % ylim([0,60]);
                % xlabel('Azimuth (Eq = 0) [°]');
                % ylabel('Elevation [°]');
                
                H = scatter3(s(:,1),s(:,2),s(:,3),20,ShR.BshF(:,1,1),'filled'); 
            else
                % OPTS = {'FaceColor','interp','edgecolor','w','edgealpha',0.2};
                OPTS = {'FaceColor','interp','edgecolor','k','linewidth',0.1,'edgealpha',0.1,...
                    'parent',ax};
                H = trisurf(T,s(:,1),s(:,2),s(:,3),ShR.BshF(:,1,1),OPTS{:});
                
                % a = solarposition.fixtoplusminus180(ShR.az-90);
                % r = sqrt(1-sind(ShR.el));
                % V = [r.*cosd(a),r.*sind(a)];
                % H = patch(ax,'Faces',T,'Vertices',V,'FaceVertexCData',mean(ShR.BshF,2:3),...
                %         'EdgeColor','none','FaceColor','interp','LineWidth',0.1);
            end
            
            hold(ax,'on'); plothalfdome(ax,1,0);
            
            K = [0.333,0.333];
            K = [eye(2),zeros(2);-diag(K)/2,diag(1+K)];
            set(ax,'position',get(ax,'position')*K);

            text(ax,0.1,0.9,'1-BshF','fontweight','bold','units','normalized');
        end
        
        function update(varargin)
        % UPDATE(T,U,P,M) - Callback function for dropdown menu changes: adjusts the way in which 
        % feature arrays F are reduced along each dimension (corresponding to menu choices), then 
        % maps the results onto the corresponding patch colors.
        
            narginchk(4,6); % args. 5 and 6 are SRC, EVENT (ignored)
            OPS = [{''} varargin(1:4)]; % let operators T,U,P,M operate on dimensions 2:5

            % Update patch colors for DwF and DshF
            for v = 1:2
                FLD = fieldnames(F{v});
                for j = 1:numel(FLD)
                    Y = reduce(F{v}.(FLD{j}),OPS);
                    if isempty(Y), continue; end
                    H{v}.(FLD{j}).CData = Y;
                end
            end
            
            % Update patch/scatter colors for BshF
            H{3}.FaceVertexCData = reduce(F{3},OPS,3:5)';
        end
        
        function Y = reduce(X,ops,dims)
        % Take 5-D array X, and apply "reducing" operations OPS(DIMS), to return a vector Y
        % that can be plotted.
        % OPS{j} = {'mean','max',.. 'std'} applies the corresponding operations along dimension j
        % OPS{j} = '42' keeps slice 42 along dimension j
        % OPS{j} = '' does nothing, the same as if size(Y,j) <= 1
            
            if nargin < 3, dims = 1:min(ndims(X),numel(ops)); end

            Y = X;
            if isempty(Y), Y = []; return; end
            for d = dims
                if size(Y,d) == 1, continue; end
                switch ops{d}
                    case '' % continue
                    case 'min', Y = min(Y,[],d);
                    case 'max', Y = max(Y,[],d);
                    case 'median', Y = median(Y,d,'omitnan');
                    case 'mode', Y = mode(Y,d);
                    case 'mean', Y = mean(Y,d,'omitnan');
                    case 'std', Y = std(Y,1,d,'omitnan');
                    otherwise, Y = arrayslice(Y,d,str2double(ops{d}));
                end
            end
        end
            
        function Y = arrayslice(X,d,p)
        % Return slice p of X along dimension d
            idx = repmat({':'},1,ndims(X));
            idx{d} = p;
            Y = X(idx{:});
        end
    end
end

%         function beamshaded(obj,varargin)
%         % R.BEAMSHADED(T,M) - returns a logical array indicating if mount(s) M are PARTIALLY-shaded
%         % (beam component) by their neighbors at simulation-step(s) T.
%         %
%         % Use R.BELOWHORIZON,
%         %
%         % The function mimics indexing over an obj.Nt×obj.Nu logical array, so that: 
%         %   R.BEAMSHADED() - returns a complete obj.Nt×obj.Nu logical array
%         %   R.BEAMSHADED(:,end) - returns an obj.Nt×1 logical vector for the last mount
%         %   R.BEAMSHADED(1:5,:) - returns steps 1:5 for all mounts.
%         %   R.BEAMSHADED(:) - returns an obj.Nt·obj.Nu vector.
%         % The only exception is linear-indexing, i.e. R.BEAMSHADED(x) for x ~= ':', that is
%         % not yet supported.
%         %
%         % See also: SHADINGRESULTS, SHADINGRESULTS.NOTDARK, SHADINGRESULTS.HSHADED
%         
%             narginchk(1,3);
%             if nargin == 2
%                 assert(varargin{1}==':','Incomplete indexing is not supported');
%                 makelist = true;
%             else, makelist = false;
%             end
%             varargin(end+1:2) = {':'}; % () and (:) to (:,:)
%             
%             idx = obj.hshmap.EIDX(varargin{1}); % zeros and valid indexes for obj.NHshaded
%             shaded = idx > 0;
%             dark = ~obj.notdark(varargin{1});
%             
%         	H(shaded,:) = obj.NHshaded(nonzeros(idx),varargin{2});  % partly shaded
%             H(dark,:) = true;                                       % all shaded when it's dark
%             H(~dark & ~shaded,:) = false;                           % none shaded
%             if makelist, H = H(:); end
%         end
%         function beamlightpolygons(obj,varargin) % {86772×1 cell}
%         end
%         function incidencefactors(obj,varargin)
%         % [Nt×Nu double]
%         end
%         function diffuseweights(obj,varargin)
%         % [Nt×Nu×4 double]
%         end
%         function diffuseshadefactors(obj,varargin)
%         % [Nt×Nu×Np×4 double]
%         end  
end % class
