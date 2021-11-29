function masks = gettrackermasks(Trackers,tolerance,IAM)
% Returns a boolean array which marks which trackers are to be considered during the shading
% analysis of each tracker in the 'analysed' index.

% There are two criteria so far, to decide wether tracker j is worth considering for the shading
% analysis of tracker i:
%	1. The maximum solid angle which j can occupy in the field of view of i must be greater than
%	   a given threshold: trackerarea/dist(i,j)^2 > 2·pi·tolerance
%	2. No other tracker of higher importance (occupying a bigger solid angle) must exist in the
%	   same direction of tracker j: acos(u(i,j)·u(i,l))^2 < 2·pi·tolerance, for any l such that
%	   dist(i,l) < dist(i,j), where u(i,k) = dist(i,k)/norm(dist(i,k));

% TODO:
%   The ful
%   - Consider angle between normals and distance-vectors, to modify area?
%   - Consider local horizon-profiles to filter out trackers behind sloping terrain.

	if nargin < 2, tolerance = getSimOption('masktolerance'); end
    
	Ntr = size(Trackers.centers,2);
	if isfield(Trackers,'analysedtrackers'), analysed = Trackers.analysedtrackers; 
	else, analysed = 1:Ntr; end
	Nu = numel(analysed);
    
	trackerarea = area(Trackers.geom.border);
    [w(1),w(2)] = rectangleproperties(Trackers.geom.border); 
    L = max(w); 
    w = min(w); % min. tracker dimension
	%distances = zeros(nat,ntr);
    
    % Calculate rotation matrix(es), for any two positions (will only be used if static)
    [F,V,~,R] = plantlayout(Trackers,[0,0],[0,90],'N2E');
    
    staticsystem = size(R,4) == 1;
    equalrotations = size(R,3) == 1;
    
    if staticsystem
    % Generate scenario...

        % Get the maximum angle [c0 = cos(th0)] between a visible point and a surface normal
        if nargin < 3, IAM = getSimOption('IAM'); end
        IAM = checkIAM(IAM);
        c0 = cosd(bisection(@(x) cosd(x).*IAM(x)-tolerance,0,90,tolerance/10));
        % TrckP = polygon3d.vf2poly(V,F);
    end
        
	masks = false(Nu,Ntr);
    waitbarh = optwaitbar(0,'Getting tracker-masks','Name','Calculating Mount-Masks...');
    for k = 1:Nu
		i = analysed(k); % i is the actual index of the analysed tracker

        d = Trackers.centers-Trackers.centers(:,i);
        r = rssq(d,1);
        r = max(w,r - L); % worst case scenario: neighbors oriented towards each other
        r(i) = Inf; % (a tracker cannot shade itself)
        
        % solid angle occupied by trackers in the field of view of tracker i > tolerance
        areacover = trackerarea./(2*pi*r.^2);
        masks(k,:) = areacover > tolerance;
        
        % minimum tracker dimension must occupy an angle larger than 32' (sun-diameter)
        masks(k,:) = masks(k,:) & w./r > 32*pi/(180*60);
        
        ic = find(masks(k,:));  % index of candidates
        d = d(:,ic)./rssq(d(:,ic),1); % reduced, normalized directions

        useful = true(size(ic));
        for j = 1:numel(ic)
            if ~useful(j), continue; end
            % seach among candidate trackers, for those in the same direction
            inline = acos(d(:,j)'*d).^2/(2*pi) < tolerance;
            inline = inline & useful;
            if areacover(ic(j)) == max(areacover(ic(inline)))
            % if the current one is the most important, remove the rest
                useful(inline) = false;
                useful(j) = true;
            else
            % otherwise remove the current one
                useful(j) = false;
            end
        end
        masks(k,ic(~useful)) = false;
        
        if staticsystem
        % Get the list of trackers which can actually be 'seen' from the analyzed tracker
        % i.e. that have any vertex more than acosd(c0) degrees in front of the tracker plane 

            idx = masks(k,:);
            if equalrotations
                n = R(:,3,1); 
            else
                n = R(:,3,i);   % surface normal of i
            end  
            x0 = V(F(i,1),:);   % a point on the surface of i
            
            x = V(F(idx,:),:); % vertices on all other mounts
            c = reshape((x-x0)*n,[],size(F,2))./r(idx)';
            idx(idx) = any(c > c0,2);

            masks(k,:) = idx;
        end

        waitbarh.update(k/Nu,sprintf('Mount %d / %d',k,Nu));
    end
end
