function varargout = pwlapprox(fun,a,b,tol,method,metric,slope,bias_K)
% [X,Y] = PWLAPPROX(FUN,A,B,TOL,METHOD,METRIC,SLOPE,BIASK) - Generates a Piece-Wise-Linear 
%   approximation to a convex/concave function FUN over interval [A,B], with an accuracy of at 
%   least TOL. The interval [A B] is split recursively at a point C = (B - A)·G + A, until the 
%   point (xc,FUN(xc)) from the original function is 'close enough' to the curve chord AB.
%
%   The partition fraction G is defined by METHOD, and the interpetation of TOL and what it means
%   for the approximation to be 'close enough' are defined by METRIC. Keyword SLOPE defines if/how
%   to use function derivatives, and BIAS_K how to attempt bias-correction. See options below.
%
%   NOTE that for non-convex/concave functions it can happen accidentally that a chord AB crosses
%   the function close enough to point C that the algorithm stops prematurely. This can be avoided
%   by splitting such functions into convex/concave segments and concatenating the results.
%
% INPUT:
%   FUN: Function handle, of the form y = @(x) FUN(x), and ideally [y,y'] = @(x) FUN(x).
%       Returning the derivative y'(x) of the function as a second argument allows error bound
%       estimation with less function evaluations, and flexibility to define partition criteria.
%       See [1] and [2], and METHOD details below.
%
%   A,B: Limits of the interval. If they're scalars they'll be interpreted as A and B, resp.
%		2-vectors are interpreted as [j,y(j)], and 3-vectors as [j,y(j),y'(j)], for j = A,B.
%		reducing the required number of function evaluations in recursive calls.
%
%   TOL: (incomplete) vector of approximation tolerances. See PARSETOLERANCE.
%
%   ERROR: character key, direction of error measurement.
%
%       'Y', let err(y) = | FUN(c) - G·y(b) + (1-G)y(a) | < tol(y) as in [1][2]
%       'X', let err(y) < | tol(x)·m |, for m = [y(b)-y(a)]/(b-a)
%       'euclidean' - let err(y) < sqrt(tol(y)² + m²·tol(x)²) 
%       'bijective' (Default) - let tol(y) < err(y) AND tol(x) < err(x)
%
%   METHOD: character key. Partition criteria, with the meaning given in [2], possibly adjusted 
%       according to the ERROR metric in use:
%
%       'midpoint' (Default) - Interval bisection (G = 1/2), doesn't require function derivatives.
%       'max' - Maximum error rule: pick the point over AB that maximizes the error between the
%           current upper- and lower-bound for FUN(x). Depending on the ERROR metric:
%               ERROR = 'Y', max. error on y: G = g0 = [y'(b) - m]/[y'(b) - y'(a)] 
%               ERROR = 'X', max. error on x: G = g0·y'(a)/m
%               otherwise, max. error perpendicular to m: G = g0·[1 + m·y'(a)]/[1+m²]
%       'chord' - Assuming that f(x) is a quadratic Bézier curve in the interval (A,B), the point
%           of maximum deviation from the chord AB (where the slope is parallel to said chord).
%           Also, the average of the 'midpoint' and 'max'-y G's: G = (1/2)(g0 + 1/2)
%       'slope' - Assuming that f(x) is a quadratic Bézier curve in the interval (A,B), the point
%           where the slope is equal to the average of the end slopes: G = (3 - 2·g0)·g0² 
%       'ensemble' - Average of all of the above partition criteria.
%
%   SLOPE: keyword {'2nd','diff','none'} For methods that require function derivatives, specify
%       whether these can be expected as the second output argument of FUN ('2nd'), or estimated
%       by differences ('diff'). Using 'none' effectively disables any method but 'bisection'.
%
%   BIAS_K: [default 0: UNSTABLE!] scalar 0 < K < 1, if K > 0, attempt to correct the bias of the 
%       resulting inner-bound PWL approximation, offseting each segment by:
% 
%       d = -K·Ey/(1-m²)·[-m;1]
%       Ey = y(b)·G - (1-G)·y(a) - y(c)
%
%     Note that d is perpendicular to the segment, and creates a vertical-offset K·Ey, where Ey 
%     is the maximum vertical error (if derivatives are available), or just the eror at the next
%     interval bisection point G. The choice of K depends on prior knowledge of the function, and
%     on the desired error metric:
% 
%       K = 0, make no correction, return points directly on the function.
%       K = 1/2 effectively centers (halves) the the maximum y-deviation (minimax optimum).
%       K = 1/3 (Default) minimizes the L2 norm of the deviations with the inner- and outer-bounds.
%       K = 1/4 is the minimax solution for a quadratic Bézier approximation.
%       ...
%
% P = PWLAPPROX(FUN,A,B,..) - For a monotonic, descending function FUN, PWLAPPROX can return an 
%   MDPWL object directly, that can be used in place of the function and its inverse:
%
%       fun(x) ~ P.val(x), fun(P.ival(y)) ~ y  within tolerance TOL 
%
% OUTPUT:
%   X,Y: Row vectors with resulting points, or...
%   P: MDPWL with the original points, P = MDPWL(X,Y,0).
%
% [1] Burkard, R.E., Hamacher, H.W., Rote, G., 1991. Sandwich approximation of univariate convex 
%   functions with an application to separable convex programming. Naval Research Logistics (NRL) 
%   (NRL) 38, 911–924. https://doi.org/10.1002/nav.3800380609
% [2] Rote, G., 1992. The convergence rate of the sandwich algorithm for approximating convex 
%   functions. Computing 48, 337–361. https://doi.org/10.1007/BF02238642
% 
% SEE ALSO: MDPWL, ONEDIODEPWLAPPROX.

    narginchk(3,8);
    % max. depth of recursive iterations, numel(x) <= 2^MaxDepth
    MaxDepth = get(0,'recursionlimit') - 1; 
    
    isnice = @(x) isfinite(x) & isreal(x);

    if nargin < 4, tol = []; end
    [~,tolx,toly] = parsetolerance(tol);
    [~,mintol] = parsetolerance(eps);

	if nargin < 5 || isempty(method), method = 'bisection'; end
    if ~ischar(method), method = 'foo'; end % crash below
    
    if nargin < 6 || isempty(metric), metric = 'euclidean'; end
    if ~ischar(metric), metric = 'foo'; end 
    
    if nargin < 7 || isempty(slope), slope = 'diff'; end
    
    % if nargin < 8 || isempty(bias_K), bias_K = 1/3; end
    if nargin < 8 || isempty(bias_K), bias_K = 0; end % bias correction is still unstable!
    assert(isscalar(bias_K) && isnice(bias_K) && bias_K >= 0 && bias_K <= 1,'Unexpected BIAS_K');
    
    slope = lower(slope);
    switch slope
    case {'none','na','off'}, hasderivatives = false;
    case {'2nd','on'}, hasderivatives = true;
    case {'diff'}
        fun = @(x) funwrap(fun,x,tolx(0)/4);
        hasderivatives = true;
    end
    
    [xa,ya,ma] = checkboundary(fun,a,hasderivatives*tolx(0)/4);
    [xb,yb,mb] = checkboundary(fun,b,hasderivatives*tolx(0)/4);
    assert(all(isnice([xa,xb,ya,yb])),'Invalid function value at edge(s) f(%s) = %s',...
    mat2str([xa,xb]),mat2str([ya,yb]));
    
    if ~all(isnice([ma,mb]))
        warning('Invalid slope values at edge(s), using bisection with no derivatives!');
        hasderivatives = false;
        method = 'bisection';
    end

    metric = lower(metric);
    switch metric
    case 'x', goodenough = @(ey,tx,ty,m) abs(ey) <= abs(tx*m);
    case 'y', goodenough = @(ey,tx,ty,m) abs(ey) <= ty;
    case 'euclidean', goodenough = @(ey,tx,ty,m) abs(ey) <= sqrt((ty^2 + (m*tx)^2)/2); 
    case 'bijective', goodenough = @(ey,tx,ty,m) abs(ey) <= min(ty,abs(m*tx));
    otherwise, error('Unexpected error metric');
    end
    
    if ~hasderivatives && ~any(strcmpi(method,{'bisection','midpoint'}))
        warning('No derivative information, using bisection method');
        method = 'bisection';
    end
    
    method = lower(method);
    switch method
    case {'bisection','midpoint'}, G = @(ma,mb,m) 0.5;
    case {'max'}
        switch metric
        case 'y', G = @gmaxEy;
        case 'x', G = @gmaxEx;
        otherwise, G = @gmaxEn; 
        end
    case {'slope'}, G = @gslope;
    case {'chord'}, G = @gchord;
    case {'ensemble'}, G = @gensemble;
    otherwise, error('Unexpected METHOD')
    end        
    assert(isnice(G(ma,mb,(yb-ya)/(xb-xa))),'Bad partition function');
    
    MaxDepthHit = false; % persistent flag for MaxDepth check
    BadGs = 0;
    try
        [x,y,y_prime] = recursivepwlapprox(xa,ya,xb,yb,ma,mb,0);
    catch ERR
        throw(ERR) % skip (big) stack of recursive calls
    end
    if BadGs > 0
       warning('pwlapprox:badGs','%d invalid partition fractions replaced by bisection',BadGs);
    end
    
    if bias_K > 0
    % Bias correction
        dx = diff(x);
        dy = diff(y);
        m = dy./dx; % segment slopes, for inner-bound approximation

        if ~hasderivatives
         % Just use the error at the midpoints (at the cost of n-1 new function evals.)
            ey_max = (y(1:end-1) + y(2:end))/2 - fun((x(1:end-1) + x(2:end))/2);
        else
            % Max. y-deviation between inner- and outer bounds
            dm = diff(y_prime);
            gg = (y_prime(2:end)-m)./dm;
            ey_max = (m-y_prime(1:end-1)).*gg.*dx;

            % ... unless slopes are too close to each other
            bad = gg < 0 | gg > 1 | abs(dm) < mintol(1) | 1./abs(dm) < mintol(1);
            if any(bad)
                idx = find(bad);    
                ey_max(idx) = (y(idx) + y(idx+1))/2 - fun((x(idx) + x(idx+1))/2);
            end
        end
        
        % Don't apply bias correction on segments that are too small
        ey_max(abs(dx) <= tolx(0) | abs(dy) <= toly(0)) = 0;

        offset_x = zeros(size(y));
        offset_y = zeros(size(y));
        dm = diff(m);

        % For each vertex j, we want to introduce an offset [dx,dy]j that corrects the bias of the
        % inner-bound-approximation, and do so in a way that offsets each segment (i,i+1) by:
        %
        %   d(i) = -K·ey_max(i)/(1-m(i)²)·[-m(i);1], for K in (0,1)          (§)
        %
        % For intermediate vertices, applying the simultaneous offset (§) on the left- and right-
        % segments results in:
        %
        %   [dx,dy]j = -K/(m(j+1)-m(j))·[1,-1;m(j+1),m(j)]·[ey_max(j);ey_max(j+1)]
        %
        offset_x(2:end-1) = bias_K*diff(ey_max)./dm;
        offset_y(2:end-1) = bias_K*(ey_max(2:end).*m(1:end-1)-ey_max(1:end-1).*m(2:end))./dm;
        
        % On the end vertices (and in cases with dm ~ 0), we just apply d(i) perpendicular to the
        % curve.
        offset_y([1,end]) = -bias_K*ey_max([1,end])./(1+m([1,end]).^2);
        offset_x([1,end]) = -m([1,end]).*offset_y([1,end]);
        
        bad = abs(dm) < mintol(1) | abs(1./dm) < mintol(1);
        if any(bad)
            idx = find(bad);
            offset_y(idx+1) = -bias_K*ey_max(idx)./(1+m(idx).^2);
            offset_x(idx+1) = -m(idx).*offset_y(idx);
        end

        if ~issorted(x + offset_x)
            % DEBUG
            t_star = x(1:end-1) + dx.*gg;
            u = y(1:end-1) + dx.*gg.*y_prime(1:end-1);
            figure(1); clf(); hold on
            plot(x,y);
            plot([x(1) t_star x(end)],[y(1) u y(end)]);
            plot(x + offset_x,y+offset_y);
            plot(x(idx),y(idx),'r.')
        end

        x = x + offset_x;
        y = y + offset_y;
    end 
    
    if nargout > 1, varargout = {x,y,y_prime};
    else
        varargout{1} = mdpwl(x,y,eps); % let mdpwl check for monotonicity
    end
    
    if MaxDepthHit
        warning('pwlapprox:MaxDepth','Recursive Iteration Limit reached, look for NaN values in [x,y]');
    end
    
    function [x,y,m] = recursivepwlapprox(xa,ya,xb,yb,ma,mb,depth)
    % Actual recursive iteration over interval (xa,ya) -> (xb,yb)
    % 'depth' keeps track of recursion depth, to check against MaxDepth
    
        mAB = (yb-ya)/(xb-xa);
        
        xc = min(abs([xa,xb])); % only used to get tolerances
        yc = min(abs([ya,yb]));
        
        if hasderivatives
        % Max y-deviation between inner- and outer bounds
            if isinf(mb)
                ey = (xb-xa)*(mAB-ma);
            else
                ey = (xb-xa)*(mAB-ma)*(mb - mAB)/(mb - ma);
            end
        else
        % Use bisection, estimate error from function value
            xc = (xa+xb)/2;
            yc = fun(xc);
            ey = (yb + ya)/2 - yc; 
        end
        
        %r = tolx(xc)/toly(yc);  % non-uniform scaling
        r = 1; 

        if goodenough(ey,tolx(xc),toly(yc),r*mAB) || abs(xb-xa) < mintol(xc) || abs(yb-ya) < mintol(yc)            
        % If the linearly interpolated value is close enough, just get out	        
            x = [xa,xb];
            y = [ya,yb];
            m = [ma,mb];
            return;
        elseif depth > MaxDepth
        % If recursion limit is reached, mark the point with NaN's
            MaxDepthHit = true;
            x = [xa,NaN,xb];
            y = [ya,NaN,yb];
            m = [ma,NaN,mb];
            return;
        end
        
            % % DEBUG
            % persistent figh
            % if depth == 0
            %     figh = {};
            %     GUIfigure('debug'); clf(); hold on;
            %     xx = linspace(xa,xb,1000);
            %     yy = fun(xx);
            %     plot(xx,yy,'b-');
            % else
            %     cellfun(@delete,figh);
            %     figh = {};
            % end
            % g0 = (mb - mAB)/(mb - ma);
            % u = [xa*(1-g0)+xb*g0,ya + (xb-xa).*g0.*ma]; % vertex of upper bound (max y-error)
            % figh{1} = plot([xa,(xa+xb)/2,xb,u(1),xa],[ya,(ya+yb)/2,yb,u(2),ya],'mo-');
            % 
            % gx = g0*ma/mAB; % max. x-error
            % figh{end+1} = plot([u(1),u(1),xa*(1-gx)+xb*gx],...
            %                    [ya*(1-g0)+yb*g0,u(2),u(2)],'go-');
            % figh{end+1} = text([u(1),xa*(1-gx)+xb*gx],[ya*(1-g0)+yb*g0,u(2)],{'My','Mx'});
            % 
            % ge = g0*(1 + mAB*ma)/(1+mAB^2); % max. Euclidean error
            % gr = g0*(r^2 + mAB*ma)/(r^2 + mAB^2);
            % figh{end+1} = plot([xa*(1-ge)+xb*ge,u(1),xa*(1-gr)+xb*gr],...
            %                    [ya*(1-ge)+yb*ge,u(2),ya*(1-gr)+yb*gr],'bo-');
            % figh{end+1} = text([xa*(1-ge)+xb*ge,xa*(1-gr)+xb*gr],...
            %                    [ya*(1-ge)+yb*ge,ya*(1-gr)+yb*gr],{'Me','Mr'});
            % 
            % 
            % % Quadratic Bézier approx.
            % t = linspace(0,1,100);
            % B = [xa;ya].*(1-t).^2 + 2*u(:).*(1-t).*t + [xb;yb].*t.^2;
            % figh{end+1} = plot(B(1,:),B(2,:),'r--');
            % 
            % t = [0.5,g0];
            % B = [xa;ya].*(1-t).^2 + 2*u(:).*(1-t).*t + [xb;yb].*t.^2;
            % figh{end+1} = plot(B(1,:),B(2,:),'ro');
            % figh{end+1} = text(B(1,:),B(2,:),{'C','S'});
            % 
            % s = [mAB,(ma + mb)/2];
            % s = 0.5*(xb - xa)*[1,1;s]./sqrt(1 + s.^2);
            % B = reshape([B-s;B+s;NaN(2)],2,6);
            % figh{end+1} = plot(B(1,:),B(2,:),'r:');
            % 
            % pause()
        
        g = G(r*ma,r*mb,r*mAB);
        if ~isfinite(g)
            g = 0.5;
            BadGs = BadGs + 1;
        elseif (xb-xa) < 2*tolx(xc)
        % Don't attempt fancy partitions on too small intervals
            g = 0.5;
        else
        % For method = 'max', metric 'x'/ 'y', and horizontal/vertical slopes,
        % G can fall too close to 0 or 1. Avoid this by setting tolerance-based boundaries
            gmax = 1-max(mintol(xb),min(tolx(xb),toly(yb)/abs(mb))/(xb-xa));
            gmin = max(mintol(xa),min(tolx(xa),toly(ya)/abs(ma))/(xb-xa));
            g = min(max(g,gmin),gmax);
        end

        xc = (xb-xa)*g+xa;
        if hasderivatives
            [yc,mc] = fun(xc);
        else
            yc = fun(xc); mc = NaN;
        end
        assert(isnice(yc),'Invalid function value: f(%f) = %f',xc,yc)

        % split the two new intervals at c, and merge the result
        [xac,yac,mac] = recursivepwlapprox(xa,ya,xc,yc,ma,mc,depth+1);
        [xcb,ycb,mcb] = recursivepwlapprox(xc,yc,xb,yb,mc,mb,depth+1);
        x = [xac, xcb(2:end)];
        y = [yac, ycb(2:end)]; % point c is skipped once, to avoid duplication
        m = [mac, mcb(2:end)];
    end
end

function [x,y,m] = checkboundary(fun,a,dx)
% Verify that the function FUN returns reasonable values at point A, and check if it returns 
% derivative information as second output.

    switch numel(a)
    case 1
        x = a;
        try
            [y,m] = fun(x);
        catch
            y = fun(x);
            m = NaN;
        end
    case 2
        x = a(1); y = a(2); 
        try
            [~,m] = fun(x);
        catch
            y = fun(x);
            m = NaN;
        end
    case 3
        x = a(1); y = a(2); m = a(3);
    otherwise
        x = NaN; y = NaN; m = NaN;
    end
    
    if ~isnan(m)
        m0 = slope(fun,x,dx);
        if abs(atand(m) - atand(m0)) > 1
           warning('Check slope function');
        end
    end
end

function [y,dy] = funwrap(fun,x,dx)
    y = fun(x);
    dy = slope(fun,x,dx);
end

function dy = slope(fun,x,dx)
% Try to get 2-sided difference approximation to fun'(x), use left- or right-sided
% difference when fun(x+dx/2) or fun(x - dx/2) are NaN or not real.
    y = [fun(x - dx/2),fun(x + dx/2)];
    if isnan(y(1)) || ~isreal(y(1))
        y(1) = fun(x);
        dx = dx/2;
    elseif isnan(y(2)) || ~isreal(y(2))
        y(2) = fun(x);
        dx = dx/2;
    end
    dy = diff(y)/dx;
end

function g = gmaxEy(ma,mb,m)  
% Partition fraction for max. y-error criterion.
    if isinf(mb), g = 1; 
    else, g = (mb - m)./(mb - ma);      
    end
end
function g = gmaxEx(ma,mb,m)  
% Partition fraction for max. x-error criterion.
    if (ma == 0), g = 0; 
    else, g = (mb - m)./(mb - ma).*ma./m;      
    end
end
function g = gmaxEn(ma,mb,m)  
% Partition fraction for max. euclidean-error criterion.
    if isinf(mb), g = (1 + m.*ma)./(1+m.^2); 
    else, g = ((mb - m)./(mb - ma)).*((1 + m.*ma)./(1+m.^2));   
    end
end
function g = gchord(ma,mb,m)  
% Partition fraction for chord criterion.
    if isinf(mb), g = 0.75; 
    else, g = 0.5*(mb - m)./(mb - ma) + 0.25;   
    end
end
function g = gslope(ma,mb,m)
% Partition fraction for slope-bisection criterion.
    g = gmaxEy(ma,mb,m);
    g = (3 - 2*g).*g.^2;
end

function g = gensemble(ma,mb,m)
% Average of all other partition criteria.
    g = gmaxEy(ma,mb,m);
    g = (g.*(1 + ma./m + (1 + m.*ma)./(1+m.^2) + (3 - 2*g).*g) + 1/2)/5;
end