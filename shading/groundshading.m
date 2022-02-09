function [shfcn,Gnd,VV,P,K] = groundshading(Trck,sunvec,Vfcn,varargin)
% [SHFCN,GND,VV,P,K] = GROUNDSHADING(TRCK,SUNVEC,VFCN) - Estimate the effect of direct shading 
%   cast over ground triangulation GND estimated for layout TRCK, given solar positions SUNVEC 
%   and plant geometry VFCN as returned by PLANTLAYOUT:
%
%         % given solar position vectors sunaz, sunel...
%         [F,VFCN,SUNVEC] = PLANTLAYOUT(TRCK,sunaz,sunel,..);
%         [SHFCN,GND,VV,P,K] = GROUNDSHADING(TRCK,SUNVEC,VFCN);
%         for t = 1:numel(sunel)
%             ...
%             % get shading of triangles GND and vertices VV at timestep t:
%             [SGND,BGND,SVV,BVV] = SHFCN(t);
%             ...
%             N = GND.faceNormal;
%             for each q (analysis_points) and R (view-point rotation)
%                
%                 W = GndVV_weights(GND,VV,P,N,q,R,PRJ,REG)
%                 REG_B = BVV'*W;
%                 REG_S = (SVV.*BVV)*'*W;
%             end
%         end
%
%  [GND,VV,P,K] are the outputs of APPROXGROUNDMESH(TRCK). GND is a rough triangulation based on
%       the TRCK. The set of vertices VV make a finer subdivision of GND (see FINEMESHPOINTS).
%       Each triangle GND(k,:) is 'sampled' at vertices VV(P(k,:),:).
%
%  Function handle SHFCN returns approximate ground "shading factors" SGND, SVV and incidence angle
%  factors BGND, BVV for each triangle GND and vertex VV, at timestep t:
%
%   [SGND,BGND,SVV,BVV] = SHFCN(t)
%
%   SGND(k), SVV(j) = 0 means no shading, i.e. incident direct irradiance over triangle k / 
%       vertex j is BHI = BNI·sind(sunel(t)).
%   SGND(k), SVV(j) = 1 means full shading, i.e. no incident direct irradiance.
%
%       NOTE! SimOptions.groundshadeangle sets a threshold (minimum angle) for which shading is 
%       actually calculated, for solar elevation angles below that threshold, all SGND and SVV 
%       are set to 1.
%
%   BGND(k) is the cosine of the solar incident angle on facet GND(k,:), i.e. the dot product of
%       the sun vector and the surface normal.
%   BVV(j) is an average of BGND(k) for all facets k that share vertex j.
%
%  The algorithm is as follows:
%       1. Project layout VFCN(t) onto GND
%       2. Let each SVV(j) = 1 if vertex VV(j,:) is shaded by any projected polygon
%       3. Estimate each SGND(k) as the weighted average of SVV(p) forall vertices VV(p,:) within
%          (or on the edge of) triangle GND(k). SGND(k) = SVV(P(k,:))·K.
%
% NOTES: (TODO) there is no self-intersection! elevations or obstacles (even if included in the 
%               mesh) will not cast shadows onto the mesh itself.
%
%        The algorithm becomes very slow at low solar elevations (long shades). At these points,
%        however, not only the relative importance of ground shading is questionable, but the
%        underlying lambertian albedo assumption most likely breaks down.
%        For these reasons there is a critical angle (SIMOPTIONS.GROUNDSHADEANGLE) below which
%        shading calculation is skipped altogether, and all SGND(k), SVV(j) are set to 1.
%
% See also: PLANTLAYOUT, PLOTTRACKERARRAY, SHADINGANALYSIS, GNDVV_WEIGHTS

    narginchk(3,Inf);
    DARK = getSimOption('groundshadeangle');
    
    validateattributes(sunvec,{'numeric'},{'real','2d','size',[NaN,3]},'','SUNVEC');
    
    parsestruct(Trck,{'centers'},'numeric','real','nonempty','size',[3,NaN]);
    parsestruct(Trck,{'geom'},'class','pvArea','scalar',@(x) isscalar(x.border));
    parsestruct(Trck,{'origin'},'numeric','real','size',[1,3]);
    Ntr = size(Trck.centers,2);

    v0 = poly2vef(Trck.geom.border,1);
    v0(:,3) = 0;
    if isfield(Trck,'axisoffset')
        validateattributes(Trck.axisoffset,{'numeric'},{'vector','finite','real','numel',3});
        v0 = v0 + Trck.axisoffset(:)';
    end
    nv = size(v0,1);
    
    try
        validateattributes(Vfcn,{'numeric','function_handle'},{'nonempty'},'','V/VFCN');
        if isnumeric(Vfcn)
            V0 = Vfcn;
            Vfcn = @(t) V0;
        else
            V0 = Vfcn(size(sunvec,1));
        end
        validateattributes(V0,{'numeric'},{'real','finite','2d','size',[Ntr*nv,3]},'','V');
    catch ERR
        error('Invalid V/VFCN: %s', ERR.message);
    end
    
    % Generate a triangulation for mount polygons. The triangulation indices will not
    % change, even when the vertices get rotated and displaced.
    Tp = delaunay(v0(:,1:2));
    Tp = repmat(Tp,Ntr,1) + repelem(size(v0,1)*(0:Ntr-1)',size(Tp,1),1);
    
    if nargin == 7 && isa(varargin{1},'triangulation')
        [Gnd,VV,P,K] = deal(varargin{4:7});
    else
        % Get rough ground mesh, and a finer set of vertices VV , for ground shading
        [Gnd,VV,P,K] = approxgroundmesh(Trck,varargin{:});
    end

    normals = Gnd.faceNormal;
    if any(normals(:,3) < 0)
       warning('Check your ground triangulation, some triangles are facing down') 
    end

    shfcn =  @arrayshading;

    function [trishade,trib,vtxshade,vtxb] = arrayshading(t)
        
        s = sunvec(t,:)';
        
        if ~all(isfinite(s))
            trishade = nan(size(Gnd,1),1);
            trib = nan(size(Gnd,1),1);
            if nargout > 2
                vtxshade = nan(size(VV,1),1);
                vtxb = nan(size(VV,1),1); 
            end
            return; 
        end
        trib = max(0,normals*s); % facet incident factors
        if nargout > 3
            vtxb = accumarray(P(:),repmat(trib,size(P,2),1),[],@mean);
        end
        
        if s(3) < sind(DARK)
            trishade = ones(size(Gnd,1),1);
            if nargout > 2
                vtxshade = ones(size(VV,1),1);
            end
            return; 
        end

        % Cast shadows onto rough Gnd mesh...
        V = Vfcn(t);
        Vp = polygon3d.projectshadow(V,s,Gnd);
        
        % ... but evaluate shading of each triangle on finer VV
        T = triangulation(Tp,Vp(:,1:2));
        vtxshade = ~isnan(T.pointLocation(double(VV(:,1:2))));
        trishade = vtxshade(P)*K';
    end
end
