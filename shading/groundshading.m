function [shfcn,Gnd,VV,P,w] = groundshading(Trck,sunvec,Vfcn,varargin)
% [SHFCN,GND,VV,P,W] = GROUNDSHADING(TRCK,SUNVEC,VFCN) - Estimate the effect of direct shading 
%   cast over ground triangulation GND estimated for layout TRCK, given solar positions SUNVEC 
%   and plant geometry VFCN as returned by PLANTLAYOUT:
%
%         % given solar position vectors sunaz, sunel...
%         [F,VFCN,SUNVEC] = PLANTLAYOUT(TRCK,sunaz,sunel,..);
%         [SHFCN,GND,VV,P,W] = GROUNDSHADING(TRCK,SUNVEC,VFCN);
%         for t = 1:numel(sunel)
%             ...
%             % get shading of triangles GND and vertices SVTX at timestep t:
%             [STRI,SVTX] = SHFCN(t);
%             ...
%         end
%
%  [GND,VV,P,W] are the outputs of APPROXGROUNDMESH(TRCK). GND is a rough triangulation based on
%       the TRCK. The set of vertices VV make a finer subdivision of GND (see FINEMESHPOINTS).
%       Each triangle GND(k,:) is 'sampled' at vertices VV(P(k,:),:).
%
%  Function handle SHFCN returns approximate ground "shading factors" [STRI,SVTX] = SHFCN(t) for
%  each triangle GND and vertex VV, at timestep t. 
%   STRI(k), SVTX(j) = 0 means no shading, i.e. incident direct irradiance over triangle k / 
%       vertex j is BHI = BNI·sind(sunel(t)).
%   STRI(k), SVTX(j) = 1 means full shading, i.e. no incident direct irradiance.
%   STRI(k), SVTX(j) < 0 is possible (!!!) when triangle GND(k,:), or one or more of the triangles
%       attached to VV(j,:) are tilted towards the sun SUNVEC(t,:). It just means the incident
%       direct irradiance on that point is higher than expected for a horizontal plane.
%
%  The algorithm is as follows:
%       1. Project layout VFCN(t) onto GND
%       2. Let each SVTX(j) = 1 if vertex VV(j,:) is shaded by any projected polygon
%       3. Estimate each STRI(k) as the weighted average of SVTX(p) forall vertices VV(p,:) within
%          (or on the edge of) triangle GND(k). STRI(k) = SVTX(P(k,:))·W.
%       4. Apply correction factors due to surface tilt: ~ max(0,n·s)/max(0.05,s(3)), where n is
%          the surface normal of each triangle, and s the sun vector at t. Correction factors for
%          vertices are taken as the average of those for all connecting triangles.
%
% NOTES: (TODO) there is no self-intersection! elevations or obstacles (even if included in the 
%               mesh) will not cast shadows onto the mesh itself.
%
%        The algorithm becomes very slow at low solar elevations (long shades). At these points,
%        however, not only the relative importance of ground shading is questionable, but the
%        underlying lambertian albedo assumption most likely breaks down.
%        For these reasons there is a critical angle (SIMOPTIONS.GROUNDSHADEANGLE) below which
%        shading calculation is skipped altogether, and all STRI(k), SVTX(j) are set to 1.
%
% See also: PLANTLAYOUT, PLOTTRACKERARRAY, SHADINGANALYSIS

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
        [Gnd,VV,P,w] = deal(varargin{4:7});
    else
        % Get rough ground mesh, and a finer set of vertices VV , for ground shading
        [Gnd,VV,P,w] = approxgroundmesh(Trck,varargin{:});
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
        trishade = vtxshade(P)*w';
    end
end
