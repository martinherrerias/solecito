
classdef mdpwl < matlab.mixin.Copyable
% Monotonic-Descending-Piece-Wise-Linear function object class. Basically an ordered set of points
% X,Y for which strict monotonicity has been enforced and thus can be used as PWL approximation to
% a bijective function: 
%
%   f(xq) ~ obj.val(xq) = interp1(obj.x,obj.y,xq,'linear')
%   f(yq) ~ obj.inval(yq) = interp1(obj.y,obj.x,yq,'linear')

properties (GetAccess = public, SetAccess = protected)
    x    % n-vector of break points
    y    % ...
end
properties (GetAccess = public, SetAccess = public)
    tol  % tolerances [dx,dy,rel] with which PWL represents underlying function
    err  % known error introduced by enforcing strict monotonicity
end
properties (Dependent = true)
    lim  % [xmin,xmax,ymin,ymax]
    n    % number of breaks
    m    % n-1 vector of segment slopes
end
methods
    function obj = mdpwl(x,y,tol,Lim)
    % PWL = MDPWL(X,Y,TOL,LIM) - Creates a Monotonic-Descending-Piece-Wise-Linear (MDPWL)
    %    interpolant for the set of points X,Y. Optionally, the curve will be simplified
    %    (i.e. 'cleaned') using a tolerance TOL, and clipped to limits LIM.
    %    As the name suggests, a strictly-monotonic, decreasing function is expected i.e.
    %    Y(k) > Y(k+1) & X(k) < X(k+1) for all k. Any deviations will try to be corrected
    %    and logged as warnings.
    %
    % X - vector of x values at break-points
    % Y - vector of y values at break-points
    % TOL - can be a vector of absolute tolerances [tx,ty], a scalar relative tolerance —   
    %     in which case tx = max(|x|·TOL,eps(x)), and similarly for ty or it can be a  
    %     3-vector [mintx,minty,reltol]. In this case, tx = max(mintx,|x|·reltol,eps(x)),  
    %     and similarly for ty. The default value is eps(0).
    %
    %     If any(TOL > 0) the points will be 'cleaned-up' before creating the MDPWL,
    %     using the Ramer-Douglas-Peucker algorithm for polygon simplification. 
    %     The algorithm can be very inefficient when there is little to clean-up, so:
    %
    %       TOL = NaN/±Inf - skip any form of input-check (faster processing, when strict
    %               monotonicity is ensured outside the function).
    %       TOL = 0 - skip Ramer-Douglas-Peucker simplification, but remove strictly 
    %               redundantpoints, and ensure strict monotonicity.
    %       TOL = eps(0) - force simplification at machine-precision.
    % 
    %     NOTE(1): Tolerances will be bound by global SimOptions.minabstol and .minreltol.
    %     For computational efficiency options are loaded only once. Use 'clear mdpwl' to
    %     reset persistent values.
    % 
    %     NOTE(2): when a relative-tolerance is provided, (i.e. scalar or 3-vector TOL), and  
    %     in order to ensure only a minimum RELATIVE-tolerace, the algorithm is applied on a
    %     quasi-logarithmic space: x' = {sgn(x)(1+ln(|x|/x0)) for |x| > x0, x/x0 otherwise},  
    %     and similarly for y0. Where x0, y0 are the minimum absolute tolerances, and the  
    %     algorithm ensures x'² + y'² < t² for a given relative-tolerance t.
    %
    % LIM - 4-vector [x_min,x_max,y_min,y_max], if provided, clips the MDPWL at the provided
    %   limits. Inf and -Inf can be used to remove bounds on one or more directions. 
    %   Lim = [], equivalent to [-Inf Inf -Inf Inf] (default) removes all bounds.
    %
    %     NOTE(3): MDPWL does not ensure that the starting and ending slopes remain strictly  
    %     unchanged (USE CAUTION WHEN EXPTRAPOLATING!).
    %
    % See also: PPVAL, INTERP1, ADDSERIES, ADDPARALLEL

        if nargin == 0, return; end % default empty object

        if nargin < 3 || isempty(tol), tol = eps(0); end     % default is machine precision
        if nargin < 4 || isempty(Lim), Lim = [-Inf Inf -Inf Inf]; end
        clipped = ~(abs(Lim)==Inf);
        
        assert(isvector(x) && isvector(y) && isreal(x) && isreal(y),'Expecting real vectors x,y');
        n = numel(x);
        assert(n > 1,'mdpwl:npts','At least 2-non-zero-points are required');

        if size(x,1) == 1, x = x'; end
        if size(y,1) == 1, y = y'; end
        assert(numel(y) == n,'mdpwl:size','Inconsistent x and y vectors');
        assert(~any(isnan(x)|isnan(y)|isinf(y)|isinf(x)),'mdpwl:NaN','NaN/Inf values in x or y');
        assert(numel(Lim) == 4,'mdpwl:lim','Limits must be a 4-vector');

        % Fast version: skip all checks if tol = Inf
        if ~any(isfinite(tol))
            obj.x = x;
            obj.y = y;
            obj.tol = NaN;
            if any(clipped), clip(obj,Lim); end
            return;
        end
        obj.tol = tol;              % remember original input
        tol = parsetolerance(tol);  % then get full tolerance vector: [ABSTOLX,ABSTOLY,RELTOL]

        % Make sure the points are sorted and monotonically decreasing
        %xy = sortrows([x,y],[1,-2]);
        xy = [x,-y];
        if any(diff(xy,1)<= 0,1), xy = unique([x,-y],'rows','sorted'); end

        obj.x = xy(:,1);
        obj.y = -xy(:,2);

        if any(clipped), clip(obj,Lim); end

        if ~all(obj.tol == 0)
        % Simplify curve using the Ramer-Douglas-Peucker algorithm on a 'tolerance-space'
            xy2tolspace(obj,tol); % overwrites xy
            xy = mdpwl.douglaspeucker([obj.x,obj.y],tol(3));
            obj.x = xy(:,1); obj.y = xy(:,2);
            xy2absspace(obj,tol); % overwrites xy
        end

        % Check that the function is in fact monotonic
        fixmonotonicity(obj,tol);

        % Resolve verticals and horizontals, to ensure strict monotonicity
        fixverticals(obj);
        fixhorizontals(obj);
    end

    function yy = val(obj,xx)
    % Y = OBJ.VAL(X) - Evaluate MDPWL function at points X.
    % See also: MDPWL.INVAL

        if isempty(xx), yy = xx; return; end

        % Check (and otherwise enforce) strict monotonicity
        if any(diff(obj.x)<=0), obj = mdpwl(obj.x,obj.y,eps(0)); end

        % get break-point intervals
        idx = mdpwl.binindex(xx,obj.x);
        
        stupidmatlab = isrow(xx);
        if stupidmatlab, xx = xx'; end  % size(x(idx)) == size(idx) EXCEPT if isrow(idx)

        % Take end slope for evaluation sites that are NaN or +-inf...
        % for obj.m(end) < 0, {NaN,Inf,-Inf} should yield {NaN,-Inf,Inf}
        nogoods = isnan(idx);
        idx(nogoods) = obj.n-1;

        % yy = zeros(size(xx),'like',xx);
        mm = (obj.y(idx+1)-obj.y(idx))./(obj.x(idx+1)-obj.x(idx));
        yy = (xx-obj.x(idx)).*mm + obj.y(idx);
        
        if stupidmatlab, yy = yy'; end
    end

    function xx = inval(obj,yy)
    % X = OBJ.INVAL(Y) - Evaluate MDPWL inverse-function at points Y, i.e. OBJ.VAL(X) = Y
    % See also: MDPWL.VAL

        if isempty(yy), xx = yy; return; end

        % Check (and otherwise enforce) strict monotonicity
        if any(diff(obj.y)>=0), obj = mdpwl(obj.x,obj.y,eps(0)); end

        % get break-point intervals
        idx = mdpwl.binindex(-yy,-obj.y);
        
        stupidmatlab = isrow(yy);
        if stupidmatlab, yy = yy'; end  % size(x(idx)) == size(idx) EXCEPT if isrow(idx)

        % Take end slope for evaluation sites that are NaN or +-inf...
        % for obj.m(end) < 0, {NaN,Inf,-Inf} should yield {NaN,-Inf,Inf}
        nogoods = isnan(idx);
        idx(nogoods) = obj.n-1;
        
        % xx = zeros(size(yy),'like',yy);
        imm = (obj.x(idx+1)-obj.x(idx))./(obj.y(idx+1)-obj.y(idx));
        xx = (yy-obj.y(idx)).*imm + obj.x(idx);
        
        if stupidmatlab, xx = xx'; end
    end

    function L = get.lim(obj)
       L = [min(obj.x),max(obj.x),min(obj.y),max(obj.y)]; 
    end
    function n = get.n(obj), n = numel(obj.x); end
    function m = get.m(obj), m = diff(obj.y)./diff(obj.x); end
    
    function yn = isvoid(obj), yn = arrayfun(@(p) numel(p.x) <= 1,obj); end

    function [Pmp,Vmp,Pmp0,Vmp0] = mpp(obj,Lim)
    % [Pmp,Vmp,Pmp0,Vmp0] = MPP(IVPWL,LIM)
    % Finds the Maximum-Power-Point, Pmp = max(v·i(v)) and Vmp = argmax(v·i(v)) for v > 0
    % on a monotonic descending current-voltage curve IVpp represented by an MDPWL object.
    % If limits are specified in LIM [v_min v_max], then the local maximum within this bounds
    % is returned in Pmp, Vmp, whereas the non-bounded MPP is returned in Pmp0, Vmp0.
    %
    % NOTE: In any case the search is limited to the first quadrant i.e. x > 0 & y > 0.
    %
    % IVPWL - MDPWL object representing a current-voltage curve, i.e. v = x, y = i
    % LIM - 2-vector [v_min v_max], Default is [0 Inf].

        % Parse input
        if nargin < 2 || isempty(Lim), Lim = [0 Inf];
        else, assert(numel(Lim)==2,'Limits must be a 2-vector');
        end
        if any(Lim < 0)
           warning('Pmp search is limited to 1Q, negative limits set to zero');
           Lim = max(Lim,0);
        end
        nolimits = Lim(1)==0 & Lim(2)==Inf;
        
        if ~isscalar(obj)
            [Pmp,Vmp,Pmp0,Vmp0] = arrayfun(@(x) mpp(x,Lim),obj);
            return;
        end

        % Start with the set of break-points in the first-quadrant
        inQ1 = obj.x >= 0 & obj.y >= 0;
        xx = obj.x(inQ1);
        yy = obj.y(inQ1);

        % Add all points for which d(xy)/dx == 0, and xi < x < xi+1 (With Local Maxima)
        % ... for an infinite line y = yi + mi(x -xi), xmp = (xi - yi/mi)/2
        % ... for xmp > 0 & ymp > 0, given mi < 0, a condition is that yi > xi·mi
        wlm = find(obj.y(1:end-1) > obj.x(1:end-1).*obj.m);
        xmp = 0.5*(obj.x(wlm)-obj.y(wlm)./obj.m(wlm));
        ininterval = xmp > obj.x(wlm) & xmp < obj.x(wlm+1);
        % first and last bounds are extended, in case MPP must be extrapolated
        if wlm(1) == 1, ininterval(1) = true; end
        if wlm(end) == numel(obj.x)-1, ininterval(end) = true; end
        xmp = xmp(ininterval);

        % make sure to include Isc, Voc, and limits (if required) - don't mind if repeated
        xlm = Lim(1:2); 
        xlm(isinf(xlm) | xlm < 0) = [];
        xlm = [xlm(:);0;obj.inval(0)];
        ylm = [obj.val(xlm(1:end-1));0];

        xx = [xx;xmp;xlm];
        yy = [yy;obj.val(xmp);ylm];

        % Find the vertex with maximum power within default boundaries (Isc to Voc)
        [Pmp,j] = max(xx.*yy);
        Vmp = xx(j);
        
        Pmp0 = Pmp; Vmp0 = Vmp;
        if nolimits, return; end

        % Check if the global MPP lies within Voltage 
        if Vmp < Lim(1) || Vmp > Lim(2)
            ininterval = xx >= Lim(1) & xx <= Lim(2);
            xx = xx(ininterval);
            yy = yy(ininterval);
            [Pmp,j] = max(xx.*yy);
            Vmp = xx(j);
        end
    end

    function clip(obj,Lim)
    % P = CLIP(P,LIM) - clip/extend the MDPWL object P to the limits LIM [x_min,x_max,y_min,y_max]. 
    %   Limits work as AND-constraints, i.e. they define a rectangular clipping window over the XY
    %   plane that the MDPWL object cannot exceed, but that it must reach. This implies that X-
    %   limits can affect Y-limits, and viceversa.
    %   Use -Inf/Inf to remove bounds on one or more directions.
    %
    % EXAMPLE:
    %   p = MDPWL([0,4],[0,-2]);      % a line with slope -1/2, p.lim = [0 4 -2 0]
    %   CLIP(p,[0 3 -1 0]);           % p.lim = [0 2 -1 0] - curve is clipped by y >= -1
    %   CLIP(p,[0 1 -1 0]);           % p.lim = [1 0 -0.5 0] - curve is clipped by x <= 1
    %   CLIP(p,[-Inf Inf -2 2]);      % p.lim = [-4 4 -2 2] - curve extends to 'touch' y = ±2
    %   CLIP(p,[0 3 -Inf Inf]);       % p.lim = [0 3 -1.5 0] - 0 <= x <= 3 are effective
    %   CLIP(p,[-Inf Inf -Inf Inf]);  % does nothing
    %
    % SEE ALSO: MDPWL

        assert(isnumeric(Lim) && isvector(Lim) && numel(Lim)==4,'Limits must be a 4-vector');
        if ~isscalar(obj)
           arrayfun(@(p) clip(p,Lim),obj);
           return
        end
        
        if obj.n == 0, return; end
        
        e = obj.n;
        clipped = [obj.x(1) < Lim(1), obj.x(e) > Lim(2), obj.y(e) < Lim(3), obj.y(1) > Lim(4)];
        extended = [obj.x(1) > Lim(1), obj.x(e) < Lim(2), obj.y(e) > Lim(3), obj.y(1) < Lim(4)];
        extended = extended & ~clipped & isfinite(Lim);
                
        if ~any(clipped | extended), return; end
        
        intrp = @(x,y,xq) (xq-x(1))./(x(2)-x(1)).*(y(2)-y(1))+y(1);
        
        % Extrapolate last segment to limits, if required
        if extended(1)
            obj.y(1) = intrp(obj.x(1:2),obj.y(1:2),Lim(1));
            obj.x(1)= Lim(1);
            clipped(4) = obj.y(1) > Lim(4);
            extended(4) = extended(4) && ~clipped(4);
        end
        if extended(2)
            obj.y(e) = intrp(obj.x(e-1:e),obj.y(e-1:e),Lim(2));
            obj.x(e)= Lim(2);
            clipped(3) = obj.y(e) < Lim(3);
            extended(3) = extended(3) && ~clipped(3);
        end
        if extended(3)
            obj.x(e) = intrp(obj.y(e-1:e),obj.x(e-1:e),Lim(3));
            obj.y(e)= Lim(3);
            clipped(2) = obj.x(e) > Lim(2);
        end
        if extended(4)
            obj.x(1) = intrp(obj.y(1:2),obj.x(1:2),Lim(4));
            obj.y(1)= Lim(4);
            clipped(1) = obj.x(1) < Lim(1);
        end
        
        if clipped(1)
        % replace the last invalid element from the start with x_min, f(x_min)
            a = find(obj.x > Lim(1),1,'first') - 1;
            obj.y(a) = intrp(obj.x(a:a+1),obj.y(a:a+1),Lim(1));
            obj.x(a)= Lim(1);    
            clipped(4) = clipped(4) && obj.y(a) > Lim(4); % is y-clipping still necessary?
        else
            a = 1;
        end
        if clipped(2)
        % replace the first invalid element from the end with x_ax, f(x_max)
            b = find(obj.x < Lim(2),1,'last') + 1;
            obj.y(b) = intrp(obj.x(b-1:b),obj.y(b-1:b),Lim(2));
            obj.x(b)= Lim(2); 
            clipped(3) = clipped(3) && obj.y(b) < Lim(3);
        else
            b = obj.n;
        end
        % ... Do the same for y bounds
        if clipped(3)
            b = find(obj.y > Lim(3),1,'last') + 1;
            obj.x(b) = intrp(obj.y(b-1:b),obj.x(b-1:b),Lim(3));
            obj.y(b)= Lim(3);
        end
        if clipped(4)
            a = find(obj.y < Lim(4),1,'first') - 1;
            obj.x(a) = intrp(obj.y(a:a+1),obj.x(a:a+1),Lim(4));
            obj.y(a)= Lim(4);
        end

        obj.x = obj.x(a:b);
        obj.y = obj.y(a:b);
    end

    function xy2tolspace(obj,tol)
    % Convert xy coordinates to 'tolerance-space', an R² space in which point-distances are 
    % obtained in units of absolute tolerance: d = (x/x0)² + (y/y0)² < 1 , for tol = [x0,y0,1]
    % or relative tolerance: d = (x'/x0)² + (y'/y0)² < t², for tol = [x0,y0,t] and t < 1
    % In this last case, it is required that: dx' = dx·{x0/x for |x| > x0, 1 otherwise},
    % which solves to x' = {sgn(x)(1+ln(|x|/x0)) for |x| > x0, x/x0 otherwise} and similarly for y

        if tol(3) <= 1
            obj.x = reltolspace(obj.x,tol(1)/tol(3));
            obj.y = reltolspace(obj.y,tol(2)/tol(3)); 
        else
            obj.x = obj.x/tol(1);
            obj.y = obj.y/tol(2);
        end

        function z = reltolspace(z,z0)
            z = z/z0;
            f = abs(z) > 1;
            z(f) = sign(z(f)).*(1+log(abs(z(f))));
        end
    end

    function xy2absspace(obj,tol)
    % Inverse of xy2tolspace(x,y,tol) 

        if tol(3) <= 1
            obj.x = abstolspace(obj.x,tol(1)/tol(3));
            obj.y = abstolspace(obj.y,tol(2)/tol(3)); 
        else
            obj.x = obj.x*tol(1);
            obj.y = obj.y*tol(2);
        end

        function z = abstolspace(z,z0)
            f = abs(z) > 1;
            z(f) = sign(z(f)).*exp(abs(z(f))-1);
            z = z*z0;
        end
    end

    function fixmonotonicity(obj,tol)
    % Use hard (i.e. x-y independent) resorting to ensure monotonicity.
    % If re-sorting implies changes beyond the preset tolerances, add them to a warning log

        if issorted(-obj.y), obj.err = 0; return; end

        ys = sort(obj.y,1,'descend');
        dy = ys - obj.y;
        dr = dy./ys;      % (y'-y)x/xy' = dy/y'
        dx = dr.*obj.x;   % (y'-y)x/y' = (dr)x
        e = nanmin(abs([dx,dy,dr])./tol,[],2);

        if any(e > 1)          
            warning('mdpwl:mono','Points not strictly sorted, RMS = %e RelTol',std(e));
            % DEBUG: ... and stop to see wtf is wrong 
            % figure(2); plot(obj.x,obj.y,'b-',xy(abs(err) > toly,1),ys(abs(err) > toly),'r.');
            % dbstack()
            % open mdpwl
            % keyboard()
            obj.err = max(abs([dx,dy,dr]),[],1);
        end
        obj.y = ys;
    end

    function fixverticals(obj)
    % Ensure strict-monotonicity by recursively replacing verticals, i.e. x(j) == x(j+1) by
    % near-verticals: x(j) - eps(x(j)) < x(j+1) + eps(x(j+1))

        tofix = diff(obj.x) <= 0;
        c = 0;
        while any(tofix)
            tofix = find(tofix);
            b = (obj.x(tofix)+obj.x(tofix+1))/2;
            obj.x(tofix) = b - eps(b);
            obj.x(tofix+1) = b + eps(b);
            tofix = diff(obj.x) <= 0;
            c = c+1;
            assert(c < 100,'fixverticals got stuck'); % DEBUG
        end
    end

    function fixhorizontals(obj)
    % Ensure strict-monotonicity by recursively replacing horizontals, i.e. y(j) == y(j+1) by
    % near-horizontals: y(j) + eps(y(j)) > y(j+1) - eps(y(j+1))

        tofix = diff(obj.y) >= 0;
        c = 0;
        while any(tofix)
            tofix = find(tofix);
            b = (obj.y(tofix)+obj.y(tofix+1))/2;
            obj.y(tofix) = b + eps(b);
            obj.y(tofix+1) = b - eps(b);
            tofix = diff(obj.y) >= 0;
            c = c+1;
            assert(c < 100,'fixverticals got stuck'); % DEBUG
        end
    end
    
    function C = interpolate(A,B,f)
    % C = INTERPOLATE(A,B,f) - Interpolate MDPWL curves A and B at fraction f using the PWL
    %   interpolation algorithm from Tsuno et al. 2009. Namely:
    %
    %   Yc = (1-f)·Ya + f·Yb'
    %   Xc = (1-f)·Xa + f·Xb'
    %   where Yb' = Ya + yb(0)-ya(0), and Xb' = xb(Yb')
    %
    % [1] Y. Tsuno, Y. Hishikawa, and K. Kurokawa, “Modeling of the I–V curves of the PV modules
    % using linear interpolation/extrapolation,” Solar Energy Materials and Solar Cells,
    % vol. 93, no. 6–7, pp. 1070–1073, Jun. 2009.

        narginchk(3,3)
        assert(isscalar(A) && isscalar(B) && isa(A,'mdpwl') && isa(B,'mdpwl'),...
            'MDPWL.INTERPOLATE requires two scalar MDPWL objects');
        assert(isscalar(f) && isreal(f) && f >= 0 && f <= 1,'f must be a scalar within [0,1]');

        % Make sure that A (the basis for Yb') is always the closest curve to C...
        if f <= 0.5
            yy = A.y + (B.val(0)-A.val(0));
            xx = B.inval(yy);
            yy = (1-f)*A.y + f*yy;
            xx = (1-f)*A.x + f*xx;
        else
        % ... if it's not, swap A <-> B and use the complement (1-f)
            yy = B.y + (A.val(0)-B.val(0));
            xx = A.inval(yy);
            yy = f*B.y + (1-f)*yy;
            xx = f*B.x + (1-f)*xx;
        end

        C = mdpwl(xx,yy,Inf);
    end
    
    function spp = scale(obj,Sx,Sy,Lim,tol)
    % S = SCALE(P,SX,SY,LIM,TOL) - scale an MDPWL object in X,Y by factors Sx, Sy, respectively.
    %   Clip and/or simplify using LIM, TOL, if required.
    %
    %   Shortcut for S = MDPWL(P.x*Sx,P.y*Sy,tol,Lim)

        narginchk(3,5);
        if nargin < 5 || isempty(tol), tol = NaN; end
        if nargin < 4 || isempty(Lim), Lim = [-Inf Inf -Inf Inf]; end
        assert(isnumeric(Sx) && isscalar(Sx) && isreal(Sx) && isfinite(Sx) &&...
               isnumeric(Sy) && isscalar(Sy) && isreal(Sy) && isfinite(Sy),...
               'Expecting scalar real factors Sx, Sy');

        spp = arrayfun(@(p) mdpwl(p.x*Sx,p.y*Sy,tol,Lim),obj);
    end

    function S = addparallel(P,Lim,tol)
    % SPP = ADDPARALLEL(IVPPS,LIM,TOL) - take two or more monotonic-descending curves (vector PP of
    %   MDPWL objects) and calculate their equivalent 'parallel' connection:
    %
    %       ys(Xs) = sum(yj(Xs)) for j = 1:M
    %
    %   Where Xs is a representative* subset of all breakpoints PP(j).x in the M original curves,  
    %   and yj(Xs) == PP(j).val(Xs) are the corresponding Y values for Xs in curve j.
    %
    %   (*) If TOL ~= NaN, the resulting MDPWL curve will be simplified to tolerance TOL; and if
    %       LIM ~= [] or [-Inf Inf -Inf Inf], it will also be clipped to limits LIM.
    %
    % INPUT:
    %   P - M-vector of MDPWL objects.
    %
    % Lim - [Vmin Vmax Imin Imax] can be specified to clip the IV curve to the required bounds.
    %		Inf and -Inf can be used to remove bounds on one or more directions. Default is no bounds.
    % Tol - is the relative tolerance (fraction of Imax) which is used as criteria to remove redundant
    %       breakpoints, default is eps(0) - actually adjusted to eps(x,y) by mkivpp 
    %
    % See also: ADDSERIES, MDPWL.MDPWL, SCALE, MPP
    
        assert(isa(P,'mdpwl') && ~isempty(P) && ~any(isvoid(P),'all'),...
            'Expecting 2+ vector of non-empty MDPWL objects');
        
        % Parse Lim = [Xmin Xmax Ymin Ymax];
        if nargin < 2 || isempty(Lim)
            Lim = [-Inf Inf -Inf Inf];
        else 
            assert(numel(Lim) == 4,'addparallel:parseLim','Limits must be a 4-vector');
        end
        if nargin < 3, tol = 0; end
        % if ~any(isfinite(tol)), tol = NaN; else, [tol,tx] = parsetolerance(tol); end

        % Reduce the list to unique curves, and keep a vector of weights (occurrences)
        [~,uidx,ic] = unique(P(:));
        W = accumarray(ic,1);
        
        M = numel(uidx);

        % For the trivial case where all curves are equal, just scale one up
        if M == 1, S = scale(P(uidx),1,W,Lim,tol); return; end
        
        % Dump vertices of unique curves into cell-array, pre-scaling Y by W
        if isvector(P) && size(P,2) > 1, W = W'; end
        XY = arrayfun(@(p,w) [p.x,w*p.y],P(uidx),W,'unif',0);
        
        % FIX! --- MEX function only works from home... COVID issues
        here = pwd();
        cd(fileparts(which('mexaddparallel')));
        
            % Call MEX function, and wrap again into MDPWL object
            [xx,yy,sz] = mexaddparallel(M,XY{:},tol,Lim);
            
        cd(here); 
        %  ---
        
        S = mdpwl(xx(1:sz)', yy(1:sz)',0);
    end
    
    function S = addseries(P,Lim,tol)
    % SPP = ADDSERIES(IVPPS,LIM,TOL) - take two or more monotonic-descending curves (vector PP of
    %   MDPWL objects) and calculate their equivalent 'series' connection:
    %
    %       xs(Ys) = sum(xj(Ys)) for j = 1:M
    %
    %   Where Ys is a representative* subset of all breakpoints PP(j).y in the M original curves,  
    %   and xj(Ys) == PP(j).inval(Ys) are the corresponding X values for Ys in curve j.
    %
    %   (*) If TOL ~= NaN, the resulting MDPWL curve will be simplified to tolerance TOL; and if
    %       LIM ~= [] or [-Inf Inf -Inf Inf], it will also be clipped to limits LIM.
    %
    % INPUT:
    %   P - M-vector of MDPWL objects.
    %
    % Lim - [Vmin Vmax Imin Imax] can be specified to clip the IV curve to the required bounds.
    %		Inf and -Inf can be used to remove bounds on one or more directions. Default is no bounds.
    % Tol - is the relative tolerance (fraction of Imax) which is used as criteria to remove redundant
    %       breakpoints, default is eps(0) - actually adjusted to eps(x,y) by mkivpp 
    %
    % See also: ADDSERIES, MDPWL.MDPWL, SCALE, MPP

        assert(isa(P,'mdpwl') && ~isempty(P) && ~any(isvoid(P),'all'),...
            'Expecting 2+ vector of non-empty MDPWL objects');
        
        % Parse Lim = [Xmin Xmax Ymin Ymax];
        if nargin < 2 || isempty(Lim)
            Lim = [-Inf Inf -Inf Inf];
        else 
            assert(numel(Lim) == 4,'addseries:parseLim','Limits must be a 4-vector');
        end
        if nargin < 3, tol = 0; end

        % Reduce the list to unique curves, and keep a vector of weights (occurrences)
        [~,uidx,ic] = unique(P(:));
        W = accumarray(ic,1);
        
        M = numel(uidx);

        % For the trivial case where all curves are equal, just scale one up
        if M == 1, S = scale(P(uidx),W,1,Lim,tol); return; end
                
        % Dump vertices of unique curves into cell-array, pre-scaling X by W
        if isvector(P) && size(P,2) > 1, W = W'; end
        XY = arrayfun(@(p,w) [w*p.x,p.y],P(uidx),W,'unif',0);
        
        % FIX! --- MEX function only works from home... COVID issues
        here = pwd();
        cd(fileparts(which('mexaddseries')));
        
            % Call MEX function, and wrap again into MDPWL object
            [xx,yy,sz] = mexaddseries(M,XY{:},tol,Lim);
            
        cd(here); 
        %  ---

        S = mdpwl(xx(1:sz)', yy(1:sz)',0);
    end

    function varargout = plot(obj,varargin)
        if ~isempty(varargin) && isa(varargin{end},'matlab.graphics.axis.Axes') && ishandle(varargin{end})
            ax = varargin{end};
            varargin(end) = [];
        else, ax = gca();
        end
        if ~ishold(), lastwill = onCleanup(@() hold('off')); end
        hold('on');
        f = ~isvoid(obj);
        h = repmat(matlab.graphics.chart.primitive.Line(),size(obj));
        h(f) = arrayfun(@(o) plot(ax,o.x,o.y,varargin{:}),obj(f));         
        if nargout > 0, varargout{1} = h; end
    end
end % methods

methods (Static = true)
    function idx = binindex(x,e)
    % IDX = binindex(X,E) - optimized equivalent to [~,IDX] = histc(X,B) for sorted x

        e(1) = -Inf;
        e(end) = Inf;
        if isscalar(x)
            idx = find(x >= e,1,'last');
        else
            [~,idx] = histc(x,e);
        end

%             a = find(x(1) >= e,1,'last');
%             if x(1)==Inf, a = a - 1; end
% 
%             if isscalar(x), idx = a; return; end
%             idx = ones(size(x));
% 
%             b = find(x(end) <= e,1,'first');
%             if x(end)==-Inf, b = 1; end
% 
%             if a == b, idx(:) = a; return; end
%             [~,idx] = histc(x,e(a:b));
%             idx = idx + a - 1;
    end

    function xy = douglaspeucker(xy,t)
    % DOUGLASPEUCKER(p,tol) - simplify the polyline represented by vector of points p (n·2 array) 
    % using the recursive Ramer-Douglas-Peucker algorithm.
    % References: http://en.wikipedia.org/wiki/Ramer%E2%80%93Douglas%E2%80%93Peucker_algorithm

        n = size(xy,1);
        if n <= 2, return; end

        %Find the point with the maximum distance
        [d2max,idx] = maxpoint2line(xy(1,:),xy(n,:),xy(2:n-1,:));
        idx = idx + 1;

        %If max distance is greater than epsilon, recursively simplify
        if sqrt(d2max) > t
            %Recursive call
            firsthalf = mdpwl.douglaspeucker(xy(1:idx,:),t);
            secondhalf = mdpwl.douglaspeucker(xy(idx:n,:),t);
            %Build the result list
            xy = [firsthalf;secondhalf(2:end,:)];
        else
            xy = xy([1,n],:); % return first and last
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
            [d2max,idx] = max(Q(:,1).*Q(:,1) + Q(:,2).*Q(:,2)); % slightly faster than sum(Q.^2,2)
        end
    end
end % static methods

end % class



