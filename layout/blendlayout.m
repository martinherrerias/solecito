function Mounts = blendlayout(Layout,Mounts)
% TRACKERS = BLENDLAYOUT(LAYOUT,MOUNTS) - Combines mount-layout information (positions and
%   orientations) with mount-model details to generate a consolidated structure fully describing
%   the physical layout of modules in a PV-plant.
%
% FUTURE: Layout and Mounts should remain separate structures, allowing the latter to be an
% array of different Tracker types (models) and the former to specify their distribution (*)
%
% Trck.name
% Trck.type {'2a','1aV','1aF','1aC','0a'} n-axis-(Fixed-tilt, Countour, Vertical), 0a = sheds
% Trck.centers - 3·Ntr array of center point coordinates
% Trck.rotation - degrees CCW from North.
% Trck.tilt - for 1aF (scalar), 0a and 1aC (Ntr-vector)
% Trck.azimuth - for 1aF (scalar), 0a and 1aC (Ntr-vector)
% Trck.slope - also 'roll angle', Ntr-vector for 0a.
% Trck.tracklimits - 2-vector for 1aX, 4-vector for 2a [azimuth,elevation]
% Trck.backtracking - boolean (for 1aF and 1aC)
% Trck.groundcoverratio - tracker width / row spacing, for back-tracking
% Trck.centerheight - (scalar) for ground-shade calculation
% Trck.axisoffset - [x,y,z] vector from center of rotation to center of tracker
% Trck.analysedtrackers - vector of Nat indices for trackers under analysis
% Trck.analysedpoints - 3·Npt array, in tracker-plane coordinate system

	narginchk(2,2);

    DEF = getSimOption('mounts');
    REQ.Mounts = {'name','type','geom','module'};
    REQ.Layout = {'name','centers','origin','rotation'};
    
    assert(isstruct(Mounts) && all(isfield(Mounts,REQ.Mounts)),...
        shortliststr(setdiff(REQ.Mounts,fieldnames(Mounts)),'Missing MOUNTS field','colon',':'));
    assert(isa(Mounts.geom,'pvArea') && Mounts.geom.depth == 4,'MOUNTS.geom must be a depth-4 pvArea');
    Np = numel(Mounts.geom.elements);

    % Check tracker-type(s) and get integer type index
    [~,trcktype] = arrayfun(@(M) ismember(lower(M.type),{'2a','1av','1af','1ac','0a'}),Mounts);
    assert(all(trcktype > 0),'blendlayout:trcktype','Unknown tracker type(s)');
    assert(isscalar(unique(trcktype)),'Multiple-mount-types are not yet supported!'); % PROVISIONAL
    trcktype = trcktype(1);
    
    % Get angle-requirements for the most-restrictive tracker type
    switch max(trcktype)
        case 1      % no angle requirements                                     % 2aV
        case 2,     REQ.Layout = [REQ.Layout,{'tilt'}];                         % 1aV
        case {3,4}, REQ.Layout = [REQ.Layout,{'tilt','azimuth'}];               % 1aF, 1aC
        case 5,     REQ.Layout = [REQ.Layout,{'tilt','azimuth','slope'}];       % 0a
    end
    
    % PROVISIONAL: merge fields into MOUNTS stucture. Allow Layout to borrow fields from Mounts, 
    % e.g. it's not clear whether 1aV-tilt is a propperty of the layout or the mount.
    for f = REQ.Layout
        if isfield(Layout,f{1}), Mounts.(f{1}) = Layout.(f{1}); end
    end
    assert(isstruct(Mounts) && all(isfield(Mounts,REQ.Layout)),...
        shortliststr(setdiff(REQ.Layout,fieldnames(Mounts)),'Missing LAYOUT field','colon',':'));
    
    % Parse centers
    noerrfun = @(x) size(x,1) == 3 && isnumeric(x) && ~any(isnan(x(:)) | isinf(x(:)));
    assert(noerrfun(Mounts.centers),'LAYOUT.centers must be a 3·N array of coordinates');
    Ntr = size(Mounts.centers,2);
    
    % Merge Info fields
    if ~isfield(Mounts,'info'), Mounts.info = {}; end
    if isfield(Layout,'info')
        Mounts.info = cat(1,Layout.info(:),Mounts.info(:));
        Mounts.info = uniquecell(Mounts.info,'stable');
    end
    
    % Set/check tracking limits, depending on tracker type
    switch trcktype
        case 1 % 2a
            Mounts = setfromdefaultifmissing(Mounts,'tracklimits',DEF.tracklimits.ver);
            noerrfun = @(x) numel(x) == 4 && x(1) <= x(2) && x(4) >= x(3) && x(4) <= 90;
            nowngfun = @(x) all((x >= [-160 80 -10 50]) & (x <= [-80 160 20 90]));
        case 2 % 1aV
            Mounts = setfromdefaultifmissing(Mounts,'tracklimits',DEF.tracklimits.ver(1:2));
            Mounts.tracklimits(3:end) = [];
            noerrfun = @(x) numel(x) == 2 && x(1) <= x(2);
            nowngfun = @(x) all((x >= [-160 80]) & (x <= [-80 160]));
        case {3,4} % 'horizontal' one-axis trackers (1aC,1aF)
            Mounts = setfromdefaultifmissing(Mounts,'tracklimits',DEF.tracklimits.hor);
            noerrfun = @(x) numel(x) == 2 && x(1) <= x(2);
            nowngfun = @(x) all((x >= [-80 20]) & (x <= [-20 80]));
        case 5 % no tracking limits
            if isfield(Mounts,'tracklimits'), Mounts = rmfield(Mounts,'tracklimits'); end
    end
    if trcktype < 5 && isscalar(Mounts.tracklimits) && Mounts.tracklimits > 0
        warning('blendlayout:scalarlims','tracklimits set to ±%0.1f',Mounts.tracklimits);
        Mounts.tracklimits = [-1,1]*Mounts.tracklimits; 
    end
    assert(trcktype == 5 || noerrfun(Mounts.tracklimits),'blendlayout:trcklims','Invalid tracking-limits');
    if trcktype ~= 5 && ~nowngfun(Mounts.tracklimits)
        warning('blendlayout:trcklims','You might want to check your tracking-limits'); 
    end

    % Find out if the system can actually move i.e. ~all(constrainedaxis)
    constrainedaxis = [trcktype == 5,trcktype ~= 1];        
    if trcktype < 5 && Mounts.tracklimits(1) == Mounts.tracklimits(2)    % for all moving things
        constrainedaxis(1) = true;
        warning('blendlayout:constr1a','1a Tracking-limits are equal, movement will be constrained');
    end
    if trcktype == 1 && Mounts.tracklimits(3)==Mounts.tracklimits(4)      % for 2a trackers
        constrainedaxis(2) = true;                  
        warning('blendlayout:constr2a','2a Tracking-limits are equal, movement will be constrained');
    end
    
    % Complete/parse LAYOUT fields
    isNtrIdx = @(x) isnumeric(x) && all(x > 0) && all(mod(x,1) == 0) && numel(x) == Ntr;

    Layout = setfromdefaultifmissing(Layout,'tidx',1:Ntr,shortliststr(1:Ntr));
    assert(isNtrIdx(Layout.tidx),'Invalid LAYOUT.tidx');
    
    Layout = setfromdefaultifmissing(Layout,'group',ones(Ntr,1),'one (single group)');
    assert(isNtrIdx(Layout.group),'Invalid LAYOUT.group');

    Layout = setfromdefaultifmissing(Layout,'model',ones(Ntr,1),'one (single model)');
    assert(isNtrIdx(Layout.model),'Invalid LAYOUT.model');
    
    % ... Build default Mount-model index from existing names
    idx = {Mounts.name}';
    lbl = strjoin(arrayfun(@(j) sprintf('%d: %s',j,idx{j}),1:numel(idx),'unif',0),', ');
    idx = LabelMap((1:numel(idx))',idx);
    Layout = setfromdefaultifmissing(Layout,'modelidx',idx,lbl);
    assert(isa(Layout.modelidx,'LabelMap') && max(Layout.model) <= prod(Layout.modelidx.esize) && ...
        numel(Layout.modelidx.pidx(Layout.model)) == Ntr,'Invalid LAYOUT.modelidx');
    
    % PROVISIONAL: merge fields into MOUNTS stucture
    Mounts.origin = Layout.origin;
    Mounts.tidx = Layout.tidx;
    Mounts.group = Layout.group;
    Mounts.model = Layout.model;
    Mounts.modelidx = Layout.modelidx;
    if isfield(Layout,'groundcoverratio'), Mounts.groundcoverratio = Layout.groundcoverratio; end
    
    % Complete/parse MOUNTS fields
    Mounts = setfromdefaultifmissing(Mounts,'axisoffset',DEF.axisoffset);
    if isscalar(Mounts.axisoffset),Mounts.axisoffset(2:3)=0; end
    noerrfun = @(x) isnumeric(x) && numel(x) == 3;
    assert(noerrfun(Mounts.axisoffset),'Invalid MOUNTS.axisoffset');
    
    Mounts = setfromdefaultifmissing(Mounts,'pidx',1:Np,shortliststr(1:Np));
    noerrfun = @(x) isnumeric(x) && isvector(x) && ~any(isnan(x) | isinf(x));
    assert(noerrfun(Mounts.pidx),'Invalid MOUNTS.pidx');

    % Set backtracking parameters for 1aC and 1aF
    if any(trcktype == [3,4]) && ~constrainedaxis(1)  
        Mounts = setfromdefaultifmissing(Mounts,'backtracking',DEF.backtracking);
        if Mounts.backtracking && ~isfield(Mounts,'groundcoverratio') 
            Mounts.groundcoverratio = guessgcr(Mounts);
            if Mounts.groundcoverratio == 0, Mounts.backtracking = 0; end
        end
        if Mounts.backtracking && (Mounts.groundcoverratio > 1 || Mounts.groundcoverratio < 0)
            error('blendlayout:gcr','groundcoverratio must lie between 0 and 1'); 
        end
    end

    if ~isfield(Mounts,'centerheight') || isempty(Mounts.centerheight)
    % Calculate center-height based on clearance from ground at critical position
  
        % Get critical sun-elevation angle
        switch trcktype
            case {3,4}, el = 90 - Mounts.tracklimits(2);  % right when backtracking kicks in
            case {2,5}, el = 0;                           % doesn't really matter
            case 1, el = max(0,Mounts.tracklimits(3));    % at lower elevation limit
        end
        
        % Calculate mean rotation matrixes at critical position (azimuth zero)
        R = mountrotations(Mounts,0,el);   % 3·3·Ntr array or 3·3 matrix
        % Keep third rows only (we only care about z)
        R = shiftdim(R(3,:,:),1); % 3·Ntr matrix 
        
        % Tracker at starting position
        p0 = polygon3d(Mounts.border);
        p0 = polytranslate(p0,Mounts.axisoffset);        
        z = poly2vef(p0)*R; % z-value of mount for all rotations
        
        Mounts.centerheight = DEF.groundclearance - min(z(:));
        warning('blendlayout:centerheight',['Pivot-height set to %0.2f m, calculated for ',...
                '(default) ground-clearance of %0.2f m'],Mounts.centerheight,DEF.groundclearance);
    end
    
    % Check row continuity for 1aC and 0a
    % if any(trcktype == [4,5]), checkangles(Mounts); end
end

function S = setfromdefaultifmissing(S,fieldname,defval,textval)
% Check structure S for field FIELDNAME, set S.FIELDNAME = DEFVAL if missing, and issue a warning
% of the form: 'S.FIELDNAME set to default TEXTVAL'

    if nargin < 3, textval = mat2str(defval); end
    if ~isfield(S,fieldname) || isempty(S.(fieldname))
        S.(fieldname) = defval;
        warning('blendlayout:missing','%s.%s set to default %s',upper(inputname(1)),fieldname,textval);
    end
end

function mgcr = guessgcr(Mounts)
% GCR = GUESSGCR(MOUNTS) - estimate plant-wide ground-cover-ratio to avoid shading, based on
% x- and z-distances of the mounts with closeby* neighbors.
% (*) The algorithm assumes near-zero azimuth, specifically, that 'important' neighbors are those
% trackers with dy < L/2, and that the shade-free angle is approx. atan2(w + dz,dx), where 
% [dx,dy,dz] is the absolute distance vector to a given neighbor, and [w, L] the tracker's width 
% and length, respectively.
%
% FUTURE: use MOUNTS.group to (optionally) assign different gcr for each row.

    [L,w] = rectangleproperties(Mounts.geom.border);
    c = Mounts.centers;
    N = size(c,2);
    
    dx = abs(c(1,:)' - c(1,:));
    dy = abs(c(2,:)' - c(2,:));
    dz = abs(c(3,:)' - c(3,:));
        
    neighbors = dy < L/2 & dx < 10*L;
    neighbors(logical(eye(N))) = false;
    if ~any(neighbors)
        warning('Could not find neighboring mounts, turning off back-tracking!');
        mgcr = 0; return;
    end
    
    gcr = nan(N);
    gcr(neighbors) = (dz(neighbors) + w)./dx(neighbors);
    gcr = nanmax(gcr,[],2);
    
    mgcr = median(gcr);
    warning(['Estimated ground-cover ratio is %0.3f-%0.3f, set to median %0.3f ',...
             '(shade-free angle of %0.1f°, and equivalent pitch %0.2f m). ',...
             'Verify value before performing shading calculations.'],...
              min(gcr),max(gcr),mgcr,asind(mgcr),w/mgcr);
end

% function checkangles(Mounts)
% % TODO: Check row continuity for 1aC and 0a
%     [L,w,offset] = rectangleproperties(Mounts.geom.border);
% 
%     N_PRE = 6;
%     THRESH = max(w,L/2);
%     S = [1/L,1/w,0]'; % dimension scale, for knn search
% 
%     c = Mounts.centers;
%     N = size(c,2);
%     
%     [R,~,Ravg] = mountrotations(Mounts,0,90);
%     S = abs(Ravg*S);
%     
%     P = polytranslate(polygon3d([-L/2,L/2],[0,0],[0,0]),Mounts.axisoffset+[offset;0]);
%     P = arrayfun(@(j) polyrotate(P,R(:,:,j)),1:size(R,3));
%     P = arrayfun(@(j) polytranslate(P(min(j,numel(P))),Mounts.centers(:,j)),1:N);
%     P = permute(reshape(poly2vef(P),2,[],3),[3 2 1]);
%     
%     p = P(:,:,1);
%     q = P(:,:,2);
%     scatter3(p(1,:),p(2,:),p(3,:)); hold on;
%     scatter3(q(1,:),q(2,:),q(3,:));
%     %%
%     neighbors = knnsearch(p,q,'K',N_PRE+1,'distance','seuclidean','scale',1./S);
%     % neighbors = full(sparse(repmat(1:N,N_PRE,1)',neighbors(:,2:end),true));
%     neighbors(:,1) = [];
%     
%     x = squeeze(R(:,1,:)); % x'-axis directions
%     q = c + x*(L/2 + offset(1));
% 
%     e = nan(N_PRE,1);
%     d_next = nan(3,N);
%     for j = N:-1:1
%         K = neighbors(j,:);
%         %d = R(:,:,j)'*(c(:,j) - c(:,ngbj)); % distance to neighbors, in x'y'z' system
%         
%         p = c(:,K) - x(:,K)*(L/2 - offset(1));
%         
%         % Unit vectors perpendicular to both x'(j) and x'(k)
%         u = cross(repmat(x(:,j),1,N_PRE),x(:,K),1);
%         u = u./rssq(u,1);
%         
%         d = rssq(p-q(:,j));
%         if any(d < THRESH)
%             e(:) = NaN;
%             for k = find(d < THRESH)
%                 e(k) = [0,0,1]*([x(:,j) -x(:,k) u(:,k)]\(q(:,j)-p(:,k)));
%             end
% 
%             [~,idx] = nanmin(abs(e(:)));
%             %if ~isnan(idx)
%                 next(j) = K(idx);
%                 d_next(:,j) = p(:,idx) - q(:,j);
%                 axis_err(:,j) = e(idx)*u(:,idx);
%             %end
%         end
% %         d = dot(q(:,j) - p,u,1);
% %         
% %         [d_next(j),idx] = min(rssq(d./S - [1;0;0],1));
% %         next(j) = K(idx);
% %         
% %         [d_prev(j),idx] = min(rssq(d./S + [1;0;0],1));
% %         prev(j) = K(idx);
%     end
%    
%     dot(d_next,x)+L;
% %     prev_ok = d_prev < THRESH & next(prev) == 1:N;
% %     next_ok = d_next < THRESH & prev(next) == 1:N;
% %     prev(~prev_ok) = NaN; next(~next_ok) = NaN;
%     %%
%     dx = atan2d(axis_err(1,:),dot(d_next,x) + L/2 + offset(1));
%     dy = atan2d(axis_err(2,:),dot(d_next,x) + L/2 + offset(1));
%     dz = atan2d(axis_err(3,:),dot(d_next,x) + L/2 + offset(1));
%     
%     rowneighbors = dx < w & dy > L/2 & dy < 1.5*L;
%     
%     neighbors = dy < L/2 & dx < 10*L;
%     neighbors(logical(eye(N))) = false;
%     if ~any(neighbors)
%         warning('Could not find neighboring mounts, turning off back-tracking!');
%         gcr = 0; return;
%     end
%     
%     gcr = nan(size(neighbors));
%     gcr(neighbors) = max((dz(neighbors) + w)./dx(neighbors));
%     [gcr,idx] = nanmax(gcr,[],2);
%     
%     d = Mounts.centers - Mounts.centers(:,idx);
%     median(atand(dx(2,:)./dx(1,:)))
%     
%     mgcr = median(gcr);
%     warning(['Estimated ground-cover ratio is %0.3f-%0.3f, set to median %0.3f ',...
%              '(shade-free angle of %0.1f°, and equivalent pitch %0.2f m). ',...
%              'Verify value before performing shading calculations'],...
%               min(gcr),max(gcr),mgcr,asind(mgcr),w/mgcr);
% end
