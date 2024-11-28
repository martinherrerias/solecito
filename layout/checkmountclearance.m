function varargout = checkmountclearance(Trackers,minclearance)
% [ZMIN, OVERLAP] = CHECKMOUNTCLEARANCE(TRACKERS) - Minimum sanity checks on physical layout, 
%   namely collision with the ground and overlap with neighbors.

    if nargin < 2, minclearance = 0; end
    
    trcktype = parselist(Trackers.type,{'2a','1aV','1aF','1aC','0a'},'mount type');
    
    % Pick a set of solar positions sunel,sunaz that could be critical for each type
    switch trcktype
        case '0a', sunel = 90; sunaz = 0;
        case {'1aF','1aC'} 
            % Make the trackers turn as far as they can go:
            % in most cases, this happens at sunel = 0° (sunset or sunrise)...
            sunel = [0,90,0]; 
            if all(isfield(Trackers,{'backtracking','groundcoverratio'})) && ...
                Trackers.backtracking && Trackers.groundcoverratio > 0
            % ... but if backtracking is used, it might happen at the shade-free-angle
                    sunel([1,3]) = asind(Trackers.groundcoverratio);
            end
            sunaz = Trackers.rotation + [-90,0,90]; % perpendicular to x' axis on 1aX trackers
        case {'2a','1aV'}
            sunaz = -180:20:180;
            sunel = zeros(size(sunaz));
    end
    
    % Generate rough ground triangular mesh:
	[v,T] = approxgroundmesh(Trackers);
    Z = v(:,3);
    Gnd = triangulation(T,v(:,1:2));
    clear v T
    
    [F,Vfcn] = plantlayout(Trackers,sunaz,sunel,'N2W');
    
    Ntr = size(F,1);
    hitsomething = false(Ntr,1);
    overlaps = zeros(Ntr,1);
    clearance = Inf(Ntr,1);
        
    for k = 1:numel(sunaz)
        if isnumeric(Vfcn), V = Vfcn; else, V = Vfcn(k); end

        % Check vertices vs interpolated ground height
        [ID,B] = Gnd.pointLocation(V(:,1:2));
        z = dot(Z(Gnd(ID,:)),B,2);
        dz = (V(:,3) - z)';
        clearance = min([clearance,dz(F)],[],2);
        
        if Ntr == 1, continue; end

        ppjh = polygon3d.vf2poly(V(:,1:2),F); % Project on XY plane
        ppjh = fixorientation(ppjh); 
        vertical = abs([ppjh.signedarea]) < eps(1);
        % Check interference with one-another

        if all(vertical), continue; end
        ppjh = pack(ppjh); % pre-pack and call clipper directly, for speed
        for j = find(~vertical)
            pji = clipper(ppjh(j),ppjh([1:j-1,j+1:end]),1,2,2);
            if ~isempty(pji)
                hitsomething(j) = true;
                pji.scale = ppjh(j).scale;
                overlaps(j) = max(overlaps(j),overlapwidth(polygon.unpack(pji)));
            end
        end
    end
    
    hitground = clearance < minclearance;
    if any(hitground)
        if ~isfield(Trackers,'tidx'),Trackers.tidx = 1:Ntr; end
        warning('checkmountclearance:gnd',[shortliststr(Trackers.tidx(hitground),'Mount',5), ...
                ' come too close to ground during operation (%0.2f m). Check centerheight ',...
                ' setting. Calculation of Diffuse shading might be seriously affected.'],...
                min(clearance));
    end
    
    if any(hitsomething)
        [overlaps,ic] = sort(overlaps,'descend');
        ic = ic(1:min(10,nnz(overlaps)));
        if ~isfield(Trackers,'tidx'),Trackers.tidx = 1:Ntr; end
        warning('checkmountclearance:trck',['%s might collide during operation, '...
            'with up to %0.0f mm overlap(s), check your layout definition!'],...
            shortliststr(Trackers.tidx(ic),'Mount',5),overlaps(1)*1000);
    end
    
    if nargout == 0, return; end
    [varargout{1:2}] = deal(clearance,overlaps);
end

function dmax = overlapwidth(p)
    for j = numel(p):-1:1
        xy = [p(j).x;p(j).y]';
        xy = xy([1:end,1],:);
        d = rssq(diff(xy),2);
        [~,idx] = max(d);
        dmax(j) = maxpoint2line(xy(idx,:),xy(idx+1,:),xy([1:idx-1,idx+2:end],:));
    end
    dmax = sqrt(max(dmax));
end

function [d2max,idx] = maxpoint2line(a,b,C)
% Calculates the maximum perpendicular distance from any point in C to the line a -> b

    Q = [C(:,1) - a(1),C(:,2) - a(2)];  % C - a
    if any(a~=b)
        u = (b-a); u = u/sqrt(u(1).^2 + u(2).^2); % unit vector in the direction a -> b
        Q = Q -(Q*u')*u;
    else
       % point-to-point distance, Q is already C - a
    end
    [d2max,idx] = max(Q(:,1).*Q(:,1) + Q(:,2).*Q(:,2));  % slightly faster than sum(Q.^2,2)
end
