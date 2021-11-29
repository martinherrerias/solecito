classdef polyprojector < matlab.mixin.Copyable
% POLYPROJECTOR - Class that holds a set of options for custom azimuthal projection of spherical 
%   polygons. Esentially each object contains a function handle (fun) for the projection, and a 
%   clipping-cone angle (acos(cX)) that determines the maximum field-of-view of the projection.  
%   There's a normalization constant r0 and a PWL approximation of the inverse function (invf).
%   The main method Q = PROJECT(P) returns a 2D polygon Q from a spherical polygon P such that:
%
%       [Q.x;Q.y] = r0·fun(P.z)·[P.x;P.y]/sqrt(P.x² + P.y²)  for  Q.z >= cX
%       P.z = invf(sqrt(Q.x² + Q.y²))
%
%   Straight lines in 3D space map to great-arcs in the unit sphere and then to potentially-complex
%   curves (e.g. elliptic, when using ortographic) in the final projection. To resolve this issue,
%   additional vertices must be added to each projected edge, to satisfy predefined approximation 
%   criteria. Properties (tol) and (angtol) control this behavior.
%
%   Option (south) allows using a mirrored projection i.e. fun(-P.z) with the same POLYPROJECTOR.
%
% PROPERTIES: split between user-accessible, flexible-input (i.e. parsing required), and 'back-end'
% transient values used for fast repetitive evaluation with the same set of settings
%
%     (GetAccess = public, SetAccess = public, AbortSet, SetObservable)
%     Event listeners for most properties* update Transient properties through LOADOBJ
%         clip*       true/false/'auto'/clipping-cone-angle (degrees)
%         normalize*  true/false/radius-of-projection for 90°, i.e. fun(0)
%         prjfun*     keyword/function-handle
%         fiam*       keyword/function-handle, integrated when prjfun is set to 'IAM'
%         tol*        relative tolerance
%         angtol      (degrees) controls edge-refinement, default sqrt(6·tol)
%         south       true/false
%
%     (Transient, GetAccess = public, SetAccess = private)
%     Updated by LOADOBJ method, based on clip, normalize, prjfun, tol and fiam
%         cX       cosine of clipping-cone angle
%         c0       cosine of 'natural' clipping-cone angle, c0 = argmax(prjfun(z)) for z in (-1:1]
%         r0       normalization constant so that r0·fun(0) == normalize
%         fun      actual (bijective & monotonic) function handle
%         invf     approx. inverse function
%
% METHODS: 
% 
%     obj = polyprojector(prjfun,varargin) - parse options and build object.
%     obj = update(obj) - triggered by public property listeners, updates transient properties
%     [Q,P] = project(obj,P,point,inside) - main function: project polygon P using current options
%     P = inverse(obj,Q) - inverse projection: get spherical polygon from flat projection Q
%
%     test(types,varargin) - generate random polygons and plot projections
%
% TODO: rewrite for face-vertex approach, to handle tessellations with redundant edges!
%       replace by existing library? e.g. https://github.com/d3/d3-geo
%
% See Also: POLYPROJECTOR.POLYPROJECTOR, POLYPROJECTOR.PROJECT, POLYGON, POLYCLIP, POLYGON3D

properties (GetAccess = public, SetAccess = public, AbortSet, SetObservable)
% Event listeners for most properties* update Transient properties through LOADOBJ
    clip        %* true/false/'auto'/clipping-cone-angle (degrees)
    normalize   %* true/false/radius-of-projection for 90°, i.e. fun(0)
    prjfun      %* keyword/function-handle
    fiam        %* keyword/function-handle, integrated when prjfun is set to 'IAM'
    tol         %* relative tolerance
    angtol      %  (degrees) controls edge-refinement, default sqrt(6·tol)
    south       %  true/false
end
properties (Transient, GetAccess = public, SetAccess = private)
% Updated by LOADOBJ method, based on clip, normalize, prjfun, tol and fiam
    cX      % cosine of clipping-cone angle (set by obj.clip)
    c0      % cosine of 'natural' clipping-cone angle, c0 = argmax(prjfun(z)) for z in (-1:1]
    r0      % normalization constant so that r0·fun(0) == normalize
    fun     % actual (bijective & monotonic) function handle
    invf    % approx. inverse function
end
properties (Hidden, GetAccess = public, SetAccess = private)
    fpp 
end
methods
    function obj = polyprojector(prjfun,varargin)
    % OBJ = POLYPROJECTOR(PRJFUN,name,value,...) - generate a POLYPROJECTOR object for projection
    % function PRJFUN, using any custom properties provided as name, value pairs, and setting the
    % rest to their default values. PRJFUN defaults to 'IAM' when empty or omitted.
    %
    % PRJFUN - function handle that returns projection radius as a function of the COSINE of the
    %   zenith angle. i.e. the +Z component of a polygon projected on the unit sphere. The function
    %   is expected to yield f(1) == 0 (i.e. +Z maps to the origin) and be monotonic increasing
    %   in [1:cX] where cX is the cosine of the clipping angle (see below). If this is not the
    %   case the function will be 'made' bijective by using a scaled azimuthal equidistant 
    %   projection for anything beyond the 'natural' clipping point c0 = argmax(prjfun(z)) 
    %
    %       fun(z) = prjfun(z)                      for z >= c0
    %       fun(z) = prjfun(c0)·acos(z)/acos(c0)    for z < c0
    %
    %   The following character-keys are also recognized as predefined functions:
    %
    %       'ort': Orthographic projection - @(z) sqrt(1-z²)
    %       'ste': Stereographic projection - @(z) sqrt(1-z²)/(1+z)
    %       'azi': Azimuthal equidistant projection - @(x) acos(z)/(pi/2)
    %       'Lam': Lambert Azimuthal equal-area projection - @(z) = sqrt(1-z)
    %
    %       'IAM': [Default] IAM-adjusted-orthographic projection. Takes Incident-Angle-Modifier
    %              function obj.fiam (see below). And integrates it to get a radiance-preserving
    %              projection, i.e. areas in the projection proportional to z·fiam(acosd(z)):
    %
    %                   @(z) integral(max(0,z·obj.fiam(acosd(z))),1,z)
    %
    %              For fast evaluation, the integral is computed once (by LOADOBJ) and saved as a
    %              griddedInterpolant by PWLAPPROX.
    %
    % POLYPROJECTOR(...,'fiam',f) - function handle or parameter structure for underlying IAM 
    %     function (when using PRJFUN = 'IAM'). See CHECKIAM for details.
    %
    % POLYPROJECTOR(...,'normalize',R0) - Scale projection so that fun(0) results in a radius R0
    %     The default is R0 = 1. Use ..,'normalize',false,.. to avoid any scaling.
    %
    % POLYPROJECTOR(...,'tol',t) - Relative tolerance (Defaults to SimOptions.RelTol). Used for PWL
    %     approximation of integrated IAM function; and for calculation of the visibility-limit TH0.
    % 
    % POLYPROJECTOR(...,'angtol',val) - straight lines in 3D space map to curves in the projection.
    %     To resolve this issue, vertices are added recursively until the projected area can be
    %     expected to be within tol·pi·(r0·fun(0))² of the 'true' projected area.
    %       
    %     The current criterium is to divide each projected arc q(A)-q(B) in segments such that
    %     acosd(q(A)·q(B))·rAB/(fun(0)·r0) < |opt.angtol|. For elliptical arcs (i.e. ortographic 
    %     projection of great-arcs), setting rAB = sqrt(r(A)·r(B)) and angtol = 360/pi·sqrt(6·tol)
    %     ensure that the projected area is within the given tolerance as long as the projected
    %     angle g = acosd(q(A)·q(B)) remains small (i.e. 1 - sin(x)/x ~ x²/6).
    %
    %     Use ..,'angtol',Inf,.. to skip this step and project only the original vertices.
    %
    % POLYPROJECTOR(...,'clip',THX) - define the angle THX (degrees) at which projections are
    %     clipped (e.g. 90° is the half dome, 180° implies no* clipping). Default ('auto' or true)
    %     is to use the largest angle TH0 = acosd(c0) for which the original projection function 
    %     is bijective, i.e. c0 = argmax(prjfun(z)) for z in (-1:1]. 'clip' = false is equivalent
    %     to THX >= 180°.
    %
    %     NOTE: Even when clip == false, and unless arctol = Inf, implicitly clipping with a cone 
    %     of sqrt(tol/2) radians is used so that the nadir [0;0;-1] is not included as part of any 
    %     arc, as it can either map to infinity or to a complete circle in the projection.
    %
    %     In any case, it is HIGHLY RECOMMENDED to pre-clip flat polygons before projection, and/or
    %     to  check that the projected polygons are not a negative (see PROJECT). 
    %
    % See Also: POLYPROJECTOR.PROJECT, POLYGON3D.CLIPWITHPLANE
    
        if nargin == 0, prjfun = 'iam'; end
        
        DEF = struct('clip','auto','normalize',true,'south',false,'prjfun',prjfun,...
          'angtol',getSimOption('angtol'),'tol',getSimOption('RelTol'),'fiam',getSimOption('IAM'));
      
        opt = getpairedoptions(varargin,DEF,'restchk');
        for f = fieldnames(opt)', obj.(f{1}) = opt.(f{1}); end
        obj = update(obj); % update cX, r0, fun

        % From here on, make listeners call LOADOBJ (MATLAB seems smart enough to avoid infinite 
        % recursive calls when setting some of these properties from inside LOADOBJ).
        addlistener(obj,'normalize','PostSet',@obj.update);
        addlistener(obj,'clip','PostSet',@obj.update);
        addlistener(obj,'prjfun','PostSet',@obj.update);
        addlistener(obj,'tol','PostSet',@obj.update);
        addlistener(obj,'fiam','PostSet',@obj.update);   
    end

    function obj = update(obj,prop,~)
    % OBJ = LOADOBJ(OBJ) -  Parse user-accessible properties to create/update transient properties:
    %   r0, cX, fun, c0, and invf.
    %
    % OBJ = LOADOBJ(OBJ,PROP) - (where PROP is a META.PROPERTY object) is meant to be called by the 
    %   event-listeners set up by the constructor upon user-accessible variables: normalize, clip,
    %   prjfun, tol, and fiam. PROP.Name is used to run only the sections of LOADOBJ that are
    %   affected by a change of PROP.
    %
    % See documentation for POLYPROJECTOR.POLYPROJECTOR.
    
        if nargin < 2, prop = []; end
        changed = @(s) isempty(prop) || any(strcmp(prop.Name,s));
        
        tol_critical = min(obj.tol,1e-6);
    
        % Defaults
        if isempty(obj.normalize) || isequal(obj.normalize,'auto'), obj.normalize = true; end
        if isempty(obj.prjfun), obj.prjfun = 'iam'; end
        
        newfunction = changed({'prjfun','tol'}) || isequal(obj.prjfun,'IAM') && changed('fiam');

        if newfunction
        % Parse/update projection function
        
            if isa(obj.prjfun,'function_handle') || isa(obj.prjfun,'griddedInterpolant')
                obj.fun = obj.prjfun;
            elseif ischar(obj.prjfun)
                switch lower(obj.prjfun(1:3))
                    case 'ste', obj.fun = @(cx) sqrt(1-cx.^2)./(1+cx);  % sin(x)/(1+cos(x))
                    case 'azi', obj.fun = @(cx) real(acosd(cx))/90;     % x/90
                    case 'lam', obj.fun = @(cx) sqrt(1-cx);             % sqrt(2)·sin(x/2)
                    case 'ort', obj.fun = @(cx) sqrt(1-cx.^2);          % sin(x)
                    case {'iam','aut'}
                        obj.fiam = checkIAM(obj.fiam);
                        obj.fpp = integrateIAM(obj.fiam,tol_critical);
                        obj.fun = @(x) obj.fpp.val(x);
                        obj.prjfun = 'IAM';
                    otherwise, obj.fun = []; % crash
                end            
            else, obj.fun = []; % crash
            end
            assert((isa(obj.fun,'function_handle') || isa(obj.fun,'griddedInterpolant')) && ...
                obj.fun(1) == 0 && obj.fun(0) > 0,'polyprojector:fun','Unrecognized/bad projection'); 
            
            % Get the 'natural' clipping point is that at which the ratio of projected areas 
            % f(cx)²/f(c90°)² == 1 within obj.tol 
            Ax_lessthan_A0 = @(x) (1-(obj.fun(x)/obj.fun(0)).^2 < 0.5*obj.tol)*1-0.5;
            obj.c0 = bisection(Ax_lessthan_A0,0,1,tol_critical);
        end

        if isnumeric(obj.normalize) && obj.normalize ~= 0
            obj.r0 = obj.normalize; 
            obj.normalize = true; 
        elseif islogical(obj.normalize)
            obj.r0 = 1;
        else
            error('polyprojector:norm','NORMALIZE must be a number (r0) or boolean value');
        end
        
        if isempty(obj.clip) || isequal(obj.clip,'auto')
        % Don't auto-clip 'well-behaved' bijective projections (e.g. lambert, azimuthal)
            % obj.clip = obj.fun(-1) <= obj.fun(0) || obj.fun(-1) > 2; % WTF?
            obj.clip = obj.fun(-1) - obj.fun(0) < obj.tol;
        end

        if changed('normalize')
            % Calculate re-normalization constant, so that projected area of hemisphere == 1
            if obj.normalize
                obj.r0 = obj.r0/obj.fun(0);
                if isnan(obj.r0) || isinf(obj.r0)
                    obj.r0 = 1;
                    warning('polyprojector:norm','Could not find normalization constant!'); 
                end
            else
                obj.r0 = 1;
            end
        end
        
        if changed('clip')
        % Calculate/update clipping-cone angle
            if obj.clip    
                if isnumeric(obj.clip) && isscalar(obj.clip) && obj.clip > 0
                    obj.cX = cosd(min(180,obj.clip));
                elseif islogical(obj.clip)
                    obj.cX = obj.c0;
                else
                    error('polyprojector:clip','CLIP option requires a logical/numeric value');
                end
            else
                obj.cX = -1;
            end
        
            if obj.cX <= -1, obj.clip = false; end
        end

        if newfunction
        % Check that the projection function is strictly monotonic, looking for argmax(r(x)) 
        % for x in [eps(1)-1,1]. Modify the function if required.
        % NOTE. eps(1)-1 avoids NaN function values e.g. for Stereographic projection.
        
            opt = optimset('TolX',tol_critical,'TolFun',0);
            [cM,~,exitflag] = fminbnd(@(x) -obj.fun(x),eps(1)-1,1,opt);
            assert(exitflag == 1,'Could not find projection maximum, check your function');
            if isfinite(obj.fun(-1)) && obj.fun(-1) > obj.fun(cM), cM = -1; end
            
            if obj.fun(obj.c0) > obj.fun(cM), cM = obj.c0; end
            % if obj.fun(cM) < obj.fun(eps(1)-1), cM = eps(1)-1; end

            % Make the cut at the natural clipping angle, if convenient
            if cM > -1 && obj.c0 > cM && 1-(obj.fun(obj.c0)/obj.fun(cM))^2 < obj.tol
                cM = obj.c0; 
            end

            if cM > eps(1)-1
            % Wrap obj.fun into a composite function to make it bijective in (-1,cM]

                if ~ischar(obj.prjfun) % skip warning for predefined functions
                    warning('polyprojector:bijective',...
                        ['A projection function that is bijective and monotonic on (-1,1] is required.'...
                        '\nUsing f(x) = r0·acos(x)/c0 for all acos(x) > acos(c0) = %0.1f°'],acosd(cM));
                end

                pfun = obj.fun; % passing obj.fun directly causes infinite recursion error
                obj.fun = @(cx) polyprojector.bijectivefcn(pfun,cM,cx);
            end
        end
        
        if newfunction || changed('clip') && obj.cX+eps(1) < obj.fpp.lim(1)
        % Calculate/update inverse function
            if isempty(obj.fpp) || any(abs(obj.fun(obj.fpp.x) - obj.fpp.y) > obj.tol)
            % Calculate/update MDPWL approximation to obj.fun
                obj.fpp = pwlapprox(obj.fun,obj.cX+eps(1),1,[1,1,1]*tol_critical);
            end
            obj.invf = @(x) obj.fpp.inval(x);
        end
        
        function p = integrateIAM(fiam,tol)
        % P = INTEGRATEIAM(FIAM) - return an interpolant (MDPWL object) for the function:
        %
        %       F(x) = integral(max(0,cos(x)*FIAM(x)),0,X)*opt.r0

            th0 = bisection(@(x) cosd(x).*fiam(x),0,90,tol);
            xx = linspace(0,th0,1/tol+1);
            % yy = cumtrapz(xx*pi/180,max(0,cosd(xx).*fiam(xx)));
            yy = sqrt(cumtrapz(xx*pi/180,sind(2*xx).*fiam(xx)));

            % f = griddedInterpolant(fliplr(cosd(xx)),fliplr(yy),'linear');
            % [cxx,yy] = pwlapprox(f,0,1,obj.tol);
            
            % % Everything beyond th0 maps to r0
            % cxx = [-1,cxx];     % cos(pi)
            % yy = [yy(1),yy];
            
            % f = griddedInterpolant(cxx,yy,'linear');
            
            % Build MDPWL p, such that p.val(cos(x)) ~ y. Note everything beyond th0 maps to r0
            p = mdpwl([cosd(xx),-1],[yy,yy(end)+eps],tol,[-1 1 -Inf Inf]);
            
            % Correct any numerical-precision offset p.val(1)
            p = mdpwl(p.x,p.y-p.val(1),0);
        end
    end
    
    function [x,r,u] = prj(obj,V)
    % [x,r,u] = PRJ(obj,V) - Project set of points V (3·N array) using polyprojector OBJ.
    %   result is the 2·N array or projected points X, their radius r and direction u.
    
        assert(isnumeric(V) && size(V,1) == 3,'Expecting 3·N array');
        V = double(V);
        
        V = V./rssq(V,1);
        r = obj.r0*obj.fun(V(3,:));
        u = V(1:2,:)./hypot(V(1,:),V(2,:));
        x = r.*u;
        x(:,r == 0) = 0;  % fix NaN values
    end

    function V = inv(obj,x)
    % V = INV(OBJ,x) - inverse of PRJ(OBJ,V), get 3D spherical coordinates from flat projection 
    %   points x (N·2 array).
        
        assert(isnumeric(x) && size(x,1) == 2,'Expecting 2·N array');
        x = double(x);
        
        r = hypot(x(1,:),x(2,:));
        cz = obj.invf(r/obj.r0);
        cz(abs(cz) > 1) = NaN;
        V = (x./r).*sqrt(1-cz.^2);
        V(3,:) = cz;
    end
    
    function [Q,P] = project(obj,P,point,inside)
    % Q = obj.project(P) - Project 3D polygon P onto the XY plane using projection function obj.fun,
    % possibly pre-clipping with a cone of angle acos(obj.cX), and scaling the result by obj.r0.
    % Since straight edges of P can map to curves in the XY plane, intermediate edges are inserted
    % and also projected (the resolution is controlled by obj.angtol).
    %
    % Even when obj.clip == false, clipping with a cone of sqrt(tol/2) radians is used to ensure
    % that the nadir [0;0;-1] is not included as part of any arc, as it can either map to infinity
    % or to a complete circle in the projection.
    %
    % Q = obj.project(P,POINT,INSIDE) - where POINT is either a 3D vector (point on space) or a 2D
    % vector (point on projection plane), and INSIDE is a boolean value, verifies and eventually
    % corrects the 'reverse-polygon' problem, by checking whether the resulting 'sign' of the
    % polygon is consistent with a POINT that is known to be INSIDE or ~INSIDE of the projection.
    %
    % FUTURE: missing overload method for face-vertex representation of sets of polygons
    %   (e.g. sky-regions, TIN surfaces) to avoid repetitive vertex projection/filtering and 
    %   edge-refining.
            
        if nargin < 3, point = []; else, point = double(point); end
        if nargin < 4, inside = true(0); else, inside = logical(inside); end
        
        if isempty(P), Q = polygon.empty; return; end
                
        arctol = sqrt(obj.tol/2); % radians

        % 0. HANDLE POLYGONS INDIVIDUALLY
        % FUTURE: an indexed-vertex-set approach would ensure that, when edges are refined (below)
        % an edge that is shared between two faces is refined only once, and its vertices projected  
        % identically for the two faces (now it might be subject to numerical errors).
        if ~isscalar(P)
            [Q,P] = arrayfun(@obj.project,P,'unif',0);
            Q = cat(1,Q{:});
            if nargout > 1, P = cat(1,P{:}); end
            return;
        end
        
        % Parse control POINT. If it falls outside the clipping cone, raise TRICKY flag, to
        % override clipping (postpone until CHECKREVERSED).
        switch numel(point)
            case 0, rp = 0;                         % no control point, not my problem
            case 2, rp = hypot(point(1),point(2));  % already projected (e.g. [0,0])
            case 3, rp = obj.r0*obj.fun(point(3));   % 3D coordinates (e.g. [0,0,1])
            otherwise, error('Expecting 2D/3D single point as POINT')
        end
            
        tricky = isnan(rp) || obj.r0*obj.fun(obj.cX) - rp <= obj.tol;
        if tricky
        % Provisionally use a well-behaved function outside cX
           fcn = obj.fun;
           lastwill = onCleanup(@() restorefunction(obj,fcn));
           obj.fun = @(cx) polyprojector.bijectivefcn(fcn,obj.cX,cx);
        end
        if numel(point) == 3, point = obj.prj(point(:)); end
        if ~all(isfinite(point))
            warning('Invalid control point, will be discarded');
            point = [];
        end

        if isvoid(P), Q = checkreversed(polygon.empty,point,inside); return; end

        A = double([P.x;P.y;P.z]);              % get list of vertices
        if obj.south, A(3,:) = -A(3,:); end
        A = A./rssq(A,1);                       % normalize, i.e. project to unit-sphere

        if obj.angtol < 180 || obj.clip || tricky
            
            % Remove any consecutive duplicate vertices (would yield 'zero-length' arcs)
            repeated = rssq(circshift(A,-1,2)-A,1) < obj.tol;
            while any(repeated)
                A(:,repeated) = [];
                % clipped(repeated) = [];
                repeated = rssq(circshift(A,-1,2)-A,1) < obj.tol;
            end
            
            if ~isempty(A)
            % 1. CLIP WITH VISIBILITY CONE 
            % Check edges for intersections with cone p·z/|p| >= cX where z is the zenith vector.
            % Insert intersection vertices for propper projection, and remove everything outside.
            % Clip also anything too close to [0;0;±1]
                [A,clipped] = clipwithcone(A);
            end
            
            if isempty(A)
            % dump invisible polygons
                Q = checkreversed(polygon.empty,point,inside); 
                P = polygon3d.empty();
                return; 
            end

            % 2. INSERT INTERMEDIATE VERTICES FOR LINE-TO-CURVE PROJECTION
            % Use recursive refinement over each great arc A(i):A(i+1)
            if obj.angtol < 180, A = refine_edges(A,~clipped); end
        end

        % 3. ACTUAL PROJECTION
        qq = obj.prj(A);
        Q = polygon(qq(1,:),qq(2,:),P.hole);
        
        % 4. CHECK/CORRECT PROJECTION SIGN        
        if ~isempty(point), Q = checkreversed(Q,point,inside); end
        
        if nargout > 1
            P.x = A(1,:); % Return refined, cliped polygon
            P.y = A(2,:);
            P.z = A(3,:)*(1-2*obj.south);
        end
        
        function restorefunction(obj,fcn), obj.fun = fcn; end

        function Q = checkreversed(Q,point,inside)
        % 4. CHECK/CORRECT PROJECTION SIGN
        
            if isempty(point), return; end
            assert(islogical(inside) && isscalar(inside),'Expecting boolean INSIDE');

            if ~isempty(Q)
                really_in = inpolygon(point(1),point(2),Q.x,Q.y);
                if Q.hole, really_in = ~really_in; end
            else, really_in = false;
            end
            
            if ~(really_in == inside)
            % Revert: add background circle and substract original projection
                % C = polygon(0:obj.angtol:360,obj.r0*obj.fun(cos(pi-arctol)),'pol');
                C = polygon(0:obj.angtol:360,obj.r0*obj.fun(obj.cX),'pol');
                Q = substractpolygons(C,Q);
            elseif tricky && obj.clip
            % (Postponed) Clipping on projection
                C = polygon(0:obj.angtol:360,obj.r0*obj.fun(obj.cX),'pol');
                Q = intersectpolygons(C,Q);
            end
        end
        
        function [A,refined] = clipwithcone(A)
        % Check edges for intersections with cone p·z/|p| >= cX where z is the zenith vector.
        % Insert intersection vertices for propper projection, and remove everything outside.
        %
        % Even when clip == false (bzw. cX <= -1), clipping is performed with double cone ±arctol 
        % to ensure that [0;0;±1] (zenith or nadir) are not included as part of any arc.
        % Nadir is problematic as it usually maps to 'anywhere on a circle' in the projection. 
        % Zenith is used to check if the projected polygon is inverted, i.e. inpolygon(zenith),
        % so we don't want it to be exactly on any edge.
        
            if ~obj.clip || tricky
            % Force clip vertices at nadir using a cone of at least arctol
            % NOTE: setting obj.cX directly avoids PostSet call to update
                cx = -cos(arctol);
            else
                cx = max(obj.cX,-cos(arctol));
            end
                           
            incone = A(3,:) >= cx;
                  
            % For each edge A -> B of P there is a plane AOB that goes through A, B, and origin O.
            % The plane AOB intersects (i.e. projects onto) the unit sphere over a great-arc AB
            % that 'goes around' plane normal w
            
            B = circshift(A,-1,2);
            w = cross(A,B);
            w = w./rssq(w,1); 
            
            % There is a possibility that some part of the arc AB crosses the cone, even when
            % both A and B lie on the same side. When this happens, insert the highest/lowest
            % point vertex, for the clipping algorithm below to detect the edge crossing:
            %
            % a) A great-circle with normal w intersects a cone with axis z and angle acosd(cX)
            %    at some point only if sqrt(1-(w·z)²) >= |cX|.
            % b) The highest/lowest points are given by x = ±uxw, where u = (zxw)/|zxw|
            % c) Points must be included in the arc AB, i.e. all(sign([a b]\x) > 0)
            % d) Don't insert new vertices if edge already crosses cone
            
            u = [0 -1 0;1 0 0;0 0 0]*w./hypot(w(1,:),w(2,:)); % (zxw)/|zxw|
            
            cxedges = xor(incone,incone([2:end,1]));  % edges that cross the cone cX (-z)
            cxarcs = ~cxedges;
            cxarcs(cxarcs) = sqrt(1-w(3,cxarcs).^2) > abs(cx);  % possible intersection

            if any(cxarcs) % || any(zxarcs)
                % Upper/lower-most arc-points, depending on whether extremes are below/above cone
                Xc = cross(u(:,cxarcs),w(:,cxarcs));
                Xc = Xc.*sign(Xc(3,:)).*(1-2*incone(cxarcs));
  
                % check if points are really within arc AB
                idx = inarc(A(:,cxarcs),B(:,cxarcs),Xc);
                cxarcs(cxarcs) = idx; 
                Xc = Xc(:,idx);
                Xc = Xc./rssq(Xc,1); % fix any numerical precision issues

                % insert (any)new vertices
                if any(cxarcs) %|| any(zxarcs)
                    [A,~] = insertvertices(A,Xc,find(cxarcs));
                    
                    incone = A(3,:) >= cx;
                    cxedges = xor(incone,incone([2:end,1]));
                end
            end
            n = size(A,2);
            
            % Two simple cases: complete polygon is gone / no clipping
            if ~any(cxedges) && ~any(incone), A = []; refined = false(0); return; end  
            if ~any(cxedges) && all(incone), refined = false(1,n); return; end  

            % Insert new vertices wherever an edge crosses the cone(s):
            if any(cxedges)
                a = find(cxedges);            
                Xc = conecrossing(cx,A(:,a),A(:,mod(a,n)+1));
                
                tooclose = any(isnan(Xc),1);
                if any(tooclose)
                    Xc = Xc(:,~tooclose);
                    a = a(~tooclose);
                end
                
                [A,new] = insertvertices(A,Xc,a);
                incone(~new) = incone;
                incone(new) = true;
                n = size(A,2);
            end

            % Now replace each cut-out region with an arc-segment of the clipping circle
            cxedges = incone-circshift(incone,-1);
            a = find(cxedges > 0);                              % first vertex of cut-out (edge)
            b = find(circshift(cxedges < 0,-a(1))) + a(1) + 1;  % last vertex of cut-out (edge)
            nc = numel(a); 

            w = zeros(3,n);    % winding direction ±[0;0;1]
            for j = 1:nc
                % ... get the flat (orthographic) projeciton of points in the cut-out region
                pxy = A(1:2,mod((a(j):b(j))-1,n)+1);
                pxy(:,all(abs(pxy) < eps(1),1)) = []; % remove point(s) exactly at zenith/nadir
                
                % ... get winding direction (±z)
                w(:,a(j)) = [0;0;sign(pxy(1,1).*pxy(2,end) - pxy(1,end).*pxy(2,1))];
                if w(3,a(j)) == 0, w(3,a(j)) = 1; end  % arbitrarily pick +z
                
                % ... if the cut-out vertices include origin, take the longer road around
                if size(pxy,2) > 2 && inpolygon(0,0,pxy(1,:),pxy(2,:))
                    w(:,a(j)) = -w(:,a(j)); 
                end
            end

            % ... remove clipped vertices
            A(:,~incone) = []; w(:,~incone) = []; cxedges(~incone) = [];
            
            % Generate vertices around clip-arcs only
            refined = cxedges > 0;
            [A,new] = refine_edges(A,refined,w);
            refined(~new) = refined;
            refined(new) = true;
        end
        
        function XX = conecrossing(cX,A,B)
        % Find the intersection(s) of line(s) A-B with cone(s) p·z/|p| = cX

            s = cross(A,B,1);             % normal* to triangle(s) OAB
            s = s./hypot(s(1,:),s(2,:));  % normalized to x²+y² = 1
            sX = sqrt(1-cX.^2);

            % Azimuth angle of cone/plane intersections (2 for each plane) - 1·N·2 array
            th = atan2d(s(2,:),s(1,:))+shiftdim([-1,1],-1).*acosd(-s(3,:).*cX./sX);

            % All cone/plane/sphere intersection points - 3·N·2 array
            XX = [sX.*cosd(th);sX.*sind(th);cX.*ones(size(th))];
            
            C = shiftdim(inarc(A,B,XX),-1);  % 1·N·2 boolean array: points really in AB?
            
            % A single (middle) point is returned when there are two intersections
            XX = sum(XX.*C,3)./sum(C,3);
            XX = XX./rssq(XX,1);
        end
        
        function C = inarc(A,B,X)
        % Take two 3·N arrays A,B and a 3·N·M array X and return an N·M boolean array C for which
        % C(i,j) = true implies that point X(:,i,j) lies within arc A(:,j), B(:,j).
        
        % X'=[A B]\X is the representation of X in terms of A, B. 
        % X' must be strictly positive for point X to be in the (shortest) arc AB
        % Using X' = [A B]\X = inv([A B]'[A B])[A B]'X, and assuming that A'A = B'B = 1
        % ... leads to the (faster!) conditions: (A'-(A·B)B')X > 0 & (A'-(A·B)B')X > 0
        
            dot = @(x,y) x(1,:,:).*y(1,:,:) + x(2,:,:).*y(2,:,:) + x(3,:,:).*y(3,:,:);

            AB = dot(A,B);           
            assert(all(abs(AB)-1 < eps(1)),'Ambiguous N·180° great-arc projection');
            C = shiftdim(dot(A - AB.*B,X) > 0 & dot(B - AB.*A,X) > 0,1); %*(1-AB.^2);

            % C = false(size(A,2),size(X,3));
            % X = permute(X,[1 3 2]);
            % for j = 1:size(A,2)
            %     % x' [a b]\x is the representation of x in terms of A, B. 
            %     % It must be strictly positive for point x to be in the (shortest) arc AB
            %     C(j,:) = all([A(:,j),B(:,j)]\X(:,:,j) > eps(1),1);
            % end
        end
        
        function [A,new] = insertvertices(A,X,p)
        % Insert new vertices X in array A after current vertices p
        
            if ~issorted(p)
               [p,idx] = sort(p);
               X = X(:,idx);
            end
            n = size(A,2); 
            nc = numel(p); 
            new = false(1,n+nc);
            new(p+(1:nc)) = true;   % positions of new vertices, in new n+nc list
            A(:,~new) = A;
            A(:,new) = X;
        end
        
        function [A,idxnew] = refine_edges(A,idxa,w,iter)
        % B = REFINE_EDGES(A) generates points over the great-arcs A(i):A(i+1) such that for
        %   each final great-arc B(j):B(j+1), either:
        %
        %       1. The angle B(j):O:B(j+1) is lower than obj.angtol, or...
        %       2. The area between the true projected arc prj(B(j):B(j+1)) and the line segment
        %          prj(B(j)):prj(B(j+1)) is estimated to be lower than obj.tol.
        %
        % [B,NIDX] = REFINE_EDGES(A,[IDX,W,ITER]) meant for recursive calls, refines only edges
        %   i:i+1 for i in find(IDX), and optionally passes precalculated plane normals W.
        %
        %   W(i) is the normal to the plane A(i)OA(i+1), W(i) = A(i)xA(i+1) /|A(i)xA(i+1)|
        %
        % Points are generated on an auxiliary coordinate system u = (A - y0), v = wxa, 
        %   y0 = (A·w)w as q(t) = u·cos(t) + v·sin(t) + y0 for t = 0...PHI.
        
            MAX_SIZE = 1e4;
            MAX_ITER = 50;
            MAX_STEP = ceil(rssq(360/obj.angtol));
        
            if nargin < 2, idxa = true(1,size(A,2)); end
            if nargin < 3
                w = cross(A,circshift(A,-1,2));
                w = w./rssq(w);
            end
            if nargin < 4, iter = 0; end

            n = size(A,2);
            
            % Divide flat projection arcs in segments g·r/r0 < |opt.angtol|
            % Each projected point q(A) will have a direction u(A) and radius r(A).
            % g is the flat angle for each arc q(A)-q(B). If the arc is approximately circular
            % i.e. r(A) ~ r(B), g·r/r0 < |angtol| = |2·sqrt(6·tol)| ensures that the projected 
            % area is within tol of the 'true' arc area.
            
            idxab = idxa | circshift(idxa,1);
            idxb = mod(find(idxa),n)+1;
            
            [q(:,idxab),r(idxab),u(:,idxab)] = obj.prj(A(:,idxab));
            
            sg = u(1,idxa).*u(2,idxb) - u(1,idxb).*u(2,idxa);   % u(i) x u(i+1) = sin(g)
            sg = sg.*sign(w(3,idxa));                           % (uxu)·(w·z)z
            g = atan2d(sg,dot(u(:,idxa),u(:,idxb)));            % signed angle difference g
            
            % Use geometric-mean radius r = sqrt(rA·rB), i.e. assume curve is ~ elliptical
            % For r < r0, g can be large while g·r/r0 remains small, and the approximation 
            % 1 - sin(x)/x ~ x²/6 no longer holds, so use r = max(sqrt(rA·rB),r0)
            rm = max(sqrt(r(idxa).*r(idxb)),obj.fun(0)*obj.r0*sqrt(obj.tol));
            ns = ceil(rm.*abs(g)/obj.angtol);
            ns(isnan(ns)) = 0;
            
            uab = q(:,idxb)-q(:,idxa);
            ds = hypot(uab(1,:),uab(2,:));
            
            ns(ds < obj.r0*sqrt(obj.tol)) = 1; 
            ns = min(ns,MAX_STEP);            % LIMIT THE REFINEMENT RATE
            idxa(idxa) = ns > 1;
            ns = ns(ns > 1);
            
            % No (further) refinement required: recursion exit
            if ~any(idxa)
                idxnew = false(1,n); 
                return; 
            end
                        
            % Points will be generated on an auxiliary coordinate system u = (A - y0), v = wxa, 
            % y0 = (A·w)w as q(t) = u·cos(t) + v·sin(t) + y0 for t = 0...PHI.

            % Get auxiliary coord. system base (only for relevant edges)
            y0 = dot(A(:,idxa),w(:,idxa)).*w(:,idxa);
            a = A(:,idxa) - y0;
            b = A(:,mod(find(idxa),n)+1) - y0;
            r = rssq(a,1);
            a = a./r;
            
            phi = dot(cross(a,b),w(:,idxa));       % (ai x aj)·w = sin(phi)
            phi = atan2d(phi,dot(a,b));
            
            v = cross(w(:,idxa),a);
            v = v./rssq(v,1);

            % ns-1 vertices will be inserted whenever ns > 1, see where to place originals
            ns = max(ns - 1,0);
            nnew = sum(ns)+n;
            
            iter = iter+1;
            if iter > MAX_ITER || nnew > MAX_SIZE
                error('Iteration/array-size limit reached')
            end
            
            idxold = zeros(1,n);
            idxold(idxa) = ns;
            idxold(2:n) = (2:n)+cumsum(idxold(1:end-1)); 
            idxold(1) = 1;
            idxnew = true(1,nnew);
            idxnew(idxold) = false;            
            A(:,~idxnew) = A;
            w(:,~idxnew) = w;

            idxold = idxold(idxa);
            
            for j = 1:nnz(idxa)
            % Only in segments where new vertices must be inserted...
            
                % Calculate actual vertices, add ns(j)-1 points linearly distributed on phi
                kth = (1:ns(j))/(ns(j)+1)*phi(j);
                newpts = r(:,j).*(a(:,j).*cosd(kth) + v(:,j).*sind(kth)) + y0(:,j);

                idxj = idxold(j)+(1:ns(j));         % where to place new vertices
                A(:,idxj) = newpts./rssq(newpts,1); % renormalization reduces numerical errors
                w(:,idxj) = repmat(w(:,idxold(j)),1,numel(idxj));
            end
            
            % Gradual, recursive mode:
            idxa(~idxnew) = idxa;
            idxa(idxnew) = true;

            [A,new] = refine_edges(A,idxa,w,iter);    % (further) refine edges
            
            idxnew(~new) = idxnew;
            idxnew(new) = true;
        end
    end
    
    function P = inverse(obj,Q)
    % P = INVERSE(OBJ,Q) - get 3D (spherical) polygon P from flat projection Q
        if ~isscalar(Q), P = arrayfun(@obj.inverse,Q); return; end
        
        V = inv(obj,[Q.x;Q.y]);
        P = polygon3d(V');
        P.hole = Q.hole;
    end
end
methods(Access = protected)
   function cpobj = copyElement(obj)
   % Called by copy(obj), propperty-listeners are instance-specific, so update after copy
        cpobj = copyElement@matlab.mixin.Copyable(obj);
        cpobj = update(cpobj);
        addlistener(cpobj,'normalize','PostSet',@cpobj.update);
        addlistener(cpobj,'clip','PostSet',@cpobj.update);
        addlistener(cpobj,'prjfun','PostSet',@cpobj.update);
        addlistener(cpobj,'tol','PostSet',@cpobj.update);
        addlistener(cpobj,'fiam','PostSet',@cpobj.update);   
    end 
end
methods (Static,Hidden = true)
    function r = bijectivefcn(pfun,c0,cx)
    % r(cx) = f(cx) for cx > c0, f(c0)·acos(cx)/acos(c0) otherwise
        r = zeros(size(cx));
        incone = cx >= c0;
        r(incone) = pfun(cx(incone));
        r(~incone) = pfun(c0).*real(acos(cx(~incone))/acos(c0));
    end
    function obj = loadobj(obj)
       obj = obj.update(); 
    end
    function test(types,varargin)
    % TEST(KEYS,...) - Compare/test projectors OBJ(k,...) for k in KEYS{:} with random polygons.
                
        if nargin == 0, types = {'Orthographic','Lambert','IAM'}; end
        if ischar(types), types = {types}; end
        assert(iscellstr(types),'Cell string of TYPES required');
 
        % Built projectors
        obj = cellfun(@(t) polyprojector(t,varargin{:}),types);
                
        noprj =  polyprojector('azim','clip',false); % used for no-clipping polygon refinement
        
        warning_reseter = naptime({'mergepolygons:arrayholes','flatten:zz'});  %#ok<NASGU>
        
        figh = gcf(); clf(); hold on;
     
        % Draw abs axis and spherical grid
        plot3([-1.3,1.3,NaN,0,0,NaN,0,0],[0,0,NaN,-1.3,1.3,NaN,0,0],[0,0,NaN,0,0,NaN,0,1.3],'k-');
        text(1.4,0,0,'X'); text(0,1.4,0,'Y'); text(0,0,1.4,'Z');
        arrayfun(@(th) polyplot(polygon3d(th,0:5:360,1,'sph'),'none',[1 1 1]*0.8),0:30:150);
        arrayfun(@(th) polyplot(polygon3d(0:5:360,th,1,'sph'),'none',[1 1 1]*0.8),-60:30:60);
        axis([-1 1 -1 1 -1 1]*1.5); axis square; axis equal;
        view(120,30);
        
        colors = hsv(2*numel(types));
        colors = colors(randi(2*numel(types),numel(types),1),:);
        colors(:,4) = 0.2;
                
        fh = {};
        while ishandle(figh)

            % Create test polygon
            F = polygon(sort(rand(8,1)*360),rand(8,1),'pol');
            
            % F = arrayfun(@(j) polygon(sort(rand(5,1)*360),rand(5,1),'pol',rand()>0.5),1:2);
            % F = mergepolygons(F);

            R = polygon3d.rotmat([rand()*180,rand()*180,rand()*180],'ZXY');
            t = rand(3,1)*2-1;
            P = polyrotate(polygon3d(F),R);
            P = polytranslate(P,t);
            
            z = [0;0;R(:,3)'*t/R(3,3)];  % the intersection of Z axis with polygon plane
            if z(3) > 0
                z = R'*(z-t);                       % ... in flat polygon coordinate system
                zin = insidepolygon(F,z(1),z(2));
            else
                zin = false;
            end
            
            [~,S] = noprj.project(P); % get refined, spherical projection

            % Delete previous projections
            cellfun(@delete,fh);
            
            fh{1} = polyplot(P,[0 0 1 0.2],'b');
            fh{2} = polyplot(S,'none','r');
            
            XY = polygon(0:5:360,1,'pol');
            for j = 1:numel(types)
                Q = obj(j).project(P,[0;0;1],zin);
                XY = substractpolygons(XY,Q);
                fh{j+3} = polyplot(Q,colors(j,:),colors(j,1:3));
            end
            fh{3} = polyplot(XY,[1 1 1 0.2]*0.8,'none');

            pause()
        end
    end
end
end
