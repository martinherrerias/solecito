function varargout = onediodepwlapprox_old(pars,Tol,Lim,g)
% [Vx,Ix] = ONEDIODEPWLAPPROX(PARS,TOL,LIM,G) - Selects a set of points Vx,Ix over the ODM curve 
%   defined by parameters PARS by recursively dividing interval LIM at a point a + G·(b-a) until 
%   a linear approximation within all intervals is reasonable within TOL. The starting interval 
%   [a,b] is provided by LIM (see below).
%
%   NOTE: If MPP, Voc, and Isc are within LIM, the algorithm makes sure they are included as
%   break points. 
%
% OBJ = ONEDIODEPWLAPPROX(PARS,TOL,LIM,G) - return an MDPWL object instead
%
%   PARS: ODM parameter set - vector [Iph Io Rs nVth Rsh] and optionally [... di2mutau Vbi Brev]
%         or structure with those field names, usually generated by TRANSLATEODM
%   LIM: 2-vector defining voltage limits for the interval. Default is [0,Voc], or...
%        4-vector of [Vmin,Vmax,Imin,Imax], equivalent to [max(Vmin,V(Imax)),min(Vmax,V(Imin))]
%        i.e. restrictions are additive, LIM defines a 'crop window'.
%   G: scalar (0..1) that defines where each interval is split. G = 0.5 (default) leads to 
%         recursive bisection.
%         Set G = (sqrt(5)-1)/2 for more elegant but less efficient golden-ratio splitting.
%   TOL: (incomplete) vector of approximation tolerances. See PARSETOLERANCE.
%
% Example: [v,i] = onediodepwlapprox([10 1e-8 0.3 1.8 400],1e-3); figure(); plot(v,i,'ro');
%          [v,i] = onediodepwlapprox([10 1e-8 0.3 1.8 400],1e-3); figure(); plot(v,i,'ro');
%
% See also: PWLAPPROX, ONEDIODEMODEL, TRANSLATEODM, MODULEINTERPOLANT, MDPWL

    NTOL = 8; % tolerance refinement factor for calls to ONEDIODEMODEL ($)

    if nargin < 4, g = 0.5; end
    if nargin < 2 || isempty(Tol), Tol = []; end
            
    % 0, Vmp, Voc are special points, that will be included whenever possible.
    [~,Vmp,Imp,Voc,Isc] = onediodeMPP(pars,0);
    
    if isempty(Tol)
        Tol = parsetolerance(Tol); % just to get Tol(3) = RelTol
        Tol = parsetolerance([Vmp,Imp,1]*Tol(3)/NTOL)*NTOL;
    end
    
    % function handle, for 'pwlapproximation' calls
    fun = @(v) onediodemodel(pars,v,Tol/NTOL); % use increased precision for onediodemodel ($)
    
    % Set interval as [0, Voc,-Inf,Inf], when omitted
    if nargin < 3, Lim = []; end
	
    switch numel(Lim)
        case 0
            LimV = [0 Voc];
            LimI = [Isc 0];
        case 2
            LimV = sort(Lim(:))';
            LimI = onediodemodel(pars,LimV,Tol/NTOL); % ($)
        case 4
            Lim(~isfinite(Lim)) = NaN;
            LimV = [Lim(1:2),onediodemodel2(pars,Lim(3:4),Tol/NTOL)];
            LimV = [max(LimV([1,4]),[],'omitnan'),min(LimV([2,3]),[],'omitnan')];
            LimI = onediodemodel(pars,LimV,Tol/NTOL);
        otherwise, LimV = [0,0]; % make it crash below
    end
    assert(LimV(2) > LimV(1),'onediodepwlapprox:lim','Bad/empty limits');
	
	% If Isc, Pmp, or Voc points are inside the interval, make sure to include them
	Vset = [0,Vmp,Voc];
	Iset = [Isc,Imp,0];
	insideinterval = Vset > LimV(1) & Vset < LimV(2);
    
    % Organize all points in a sorted array of points
    VIset(:,1) = [Vset(insideinterval),LimV]';
    VIset(:,2) = [Iset(insideinterval),LimI]';
    VIset = unique(VIset,'rows')'; % ends up as a 2xn array [Vx;Ix], where 2 <= n <= 4
	
    %clf(); plot(VIset(1,:),VIset(2,:),'bo'); hold on; %DEBUG!
    
	% Split all intervals, and concatenate results
	vv = cell(4,1);
	ii = cell(4,1);
	for j = 1:size(VIset,2)-1
		[vv{j},ii{j}] = pwlapprox(fun,VIset(:,j),VIset(:,j+1),Tol,g);
	end
	Vx = [vv{1},vv{2},vv{3},vv{4}];
	Ix = [ii{1},ii{2},ii{3},ii{4}]; % unused intervals will be empty anyway
 
    % GUIfigure('debug'); clf();
    % plot(Vx,Ix)
    % hold on
    % plot(Vx,Ix,'r.')
    % foo = 0;

    switch nargout
        case 2, varargout = {Vx,Ix};
        case {0,1}, varargout = {mdpwl(Vx,Ix,eps(0))};
        otherwise, error('Too many output arguments')
    end
end
