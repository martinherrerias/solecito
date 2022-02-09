function VW = GndVV_weights(Gnd,GndVV,GndP,Gnd_normals,Pref0,Rpov,IAMprj,Pj_Gnd)
% W = GNDVV_WEIGHTS(GND,VV,P,N,Q,R,PRJ,REG) - Estimate integration weights for the TIN of ground
%   vertices VV,P, as returned by GROUNDSHADING, for a given view-point Q, R and a projection
%   PRJ. The result W is a size(VV,1) x numel(REG) sparse matrix with non-zero weights W(j,k) for
%   every vertex VV(j) that is visible from Q and belongs to projected region REG{k}
%
%         % given solar position vectors sunaz, sunel...
%         [F,VFCN,SUNVEC] = PLANTLAYOUT(TRCK,sunaz,sunel,..);
%         [SHFCN,GND,VV,P] = GROUNDSHADING(TRCK,SUNVEC,VFCN);
%         for t = 1:numel(sunel)
%             ...
%             % get shading of triangles GND and vertices VV at timestep t:
%             [SGND,BGND,SVV,BVV] = SHFCN(t);
%             ...
%             N = GND.faceNormal;
%             for each q (analysis_points) and R (view-point rotation)
%                
%                 W = GNDVV_WEIGHTS(GND,VV,P,N,q,R,PRJ,REG)
%                 REG_B = BVV'*W;
%                 REG_S = (SVV.*BVV)*'*W;
%             end
%         end
%
%   TODO: after clipping the triangulation, it should be propperly projected and merged with the
%   projected regions, to account for (and return) local horizons.
%
% See also: SHADINGANALYSIS, GROUNDSHADING, PLANTLAYOUT

    % Rotate rough ground mesh, and reduce to visible facets
    GndV_POA = (Gnd.Points - Pref0')*Rpov';
    GndT_POA = Gnd.ConnectivityList;

    visible_pts = GndV_POA(:,3) > IAMprj.cX;

    if ~any(visible_pts)
        VW = spalloc(size(GndVV,1),1,0);
        return;
    end
    
    visible_tri = any(visible_pts(GndT_POA),2);
    point_away = visible_tri;
    point_away(point_away) = ...
        dot(Pref0'-Gnd.Points(GndT_POA(point_away,1),:),Gnd_normals(point_away,:),2) <= 0;
    visible_tri = visible_tri & ~point_away;

    % Rotate fine mesh points (only for visible facets)
    [iv,~,GndP_POA] = unique(GndP(visible_tri,:));
    GndP_POA = reshape(GndP_POA,[],size(GndP,2));
    GndVV_POA = (GndVV(iv,:) - Pref0')*Rpov';                

    % Project mesh points, and assign vertex weights based on the projected area of 
    % their connected triangles (weight zero for points behind visibility cone).
    [VVprj,r] = IAMprj.prj(GndVV_POA');
    % W = approxprjarea(VVprj',GndP_POA(:,1:3));
    W = approxprjarea(GndVV_POA,GndP_POA(:,1:3),IAMprj.cX);

    vtxW = accumarray(GndP_POA(:),repmat(W,size(GndP_POA,2),1),[],@mean);
    vtxW = vtxW.*sqrt(1-r'.^2);
    vtxW(GndVV_POA(:,3) <= IAMprj.cX) = 0;
    
    visible_gnd = ~cellfun(@isempty,Pj_Gnd);
    VW = repmat({spalloc(size(GndVV,1),1,0)},1,numel(visible_gnd));

    for j = find(visible_gnd)
        inj = polygon.inpackedpolygon(Pj_Gnd{j},VVprj(1,:),VVprj(2,:));
        if any(inj)
            w = vtxW(inj);
            w = w/sum(w);
            % vtxW(inj) = w;
            VW{j} = sparse(iv(inj),1,w,size(GndVV,1),1);
        end
    end
    VW = cat(2,VW{:});
end

function A = approxprjarea(V,T,h)

    if nargin < 3, h = 0; end
   
    [iu,~,T] = unique(T);
    V = V(iu,:);
    T = reshape(T,[],3);
    
    Vp = V./rssq(V,2);
    below = Vp(:,3) <= h;
    type = sum(below(T),2);
    
    % find all edges that cross XY plane
    E = unique(sort([T(:,1:2);T(:,2:3);T(:,[1,3])],2),'rows','sorted');
    E = E(sum(below(E),2) == 1,:);
    
    % add their intersection point to the list of points
    f = V(E(:,1),3)./(V(E(:,1),3) - V(E(:,2),3));
    q = V(E(:,2),:).*f + V(E(:,1),:).*(1-f);
    iq = size(V,1)+(1:size(q,1))';
    
    V = [Vp;q./rssq(q,2)];
    r = hypot(V(:,1),V(:,2));
    
    A = zeros(size(T,1),1);
    
    for typ = 0:2
        filter = (type == typ);
        if ~any(filter), continue; end
        
        Tft = T(filter,:)';
        
        switch typ
        case 0
        % all points above horizon

        case 1
        % one vertex below horizon: replace by rectangle (two edge intersections)
            
            tochange = Tft(below(Tft));
            for l = 1:3
                isatend = tochange == Tft(3,:)';
                Tft(:,~isatend) = circshift(Tft(:,~isatend),1,1);
            end
            R = sort(Tft([3,1],:),1)';
            idx = knnsearch(E,R);
            Tft(4,:) = iq(idx);
            
            R = sort(Tft(2:3,:),1)';
            idx = knnsearch(E,R);
            Tft(3,:) = iq(idx);

        case 2
        % two vertices below horizon: replace vertices by edge intersections
        
            tokeep = Tft(~below(Tft));
            for l = 1:3
                isatend = tokeep == Tft(3,:)';
                Tft(:,~isatend) = circshift(Tft(:,~isatend),1,1);
            end
            R = sort(Tft([3,1],:),1)';
            idx = knnsearch(E,R);
            Tft(1,:) = iq(idx);
            
            R = sort(Tft(2:3,:),1)';
            idx = knnsearch(E,R);
            Tft(2,:) = iq(idx);

        % case 3 
        %   all points below horizon
        end
        
        for l = 1:size(Tft,1)
            next = mod(l,size(Tft,1))+1;
            A(filter) = A(filter) + segment(Tft(l,:)',Tft(next,:)');
        end
    end
    
    A = abs(A);

    function a = segment(jj,kk)
        rjk = r(jj).*r(kk);
        J = V(jj,:);
        K = V(kk,:);
        c = dot(J,K,2);
        s = J(:,1).*K(:,2)-J(:,2).*K(:,1);
        a = rjk.*atan2(s,c);
    end
end