function [Q,W] = sunposmesh(varargin)
% [V,W,PRJ] = SUNPOSMESH(S) - Generate a triangular mesh of points V over the unit sphere, that
%   fully contains the set of unit vectors S, and can be used to interpolate the values of a
%   scalar-valued function f(S) as f(S) ~ W·f(V) from a weight matrix W and the reduced set V.
%
%   Weights are calculated by W = INTERPMATRIX(S,V,'-sph').
%
% [V,W,PRJ] = SUNPOSMESH(AZ,EL) - Use azimuth/elevation vectors AZ, EL (degrees) instead of S.
%
% [V,W,PRJ] = SUNPOSMESH(..,'name',value): See ['option',default] pairs below.
%
%     'meshsize',SimOptions.angtol·pi/180 - (scalar) max. projected mesh edge length ; or
%       mesh-size function handle f = @(xy) that takes projected points xy and returns max edge 
%       lengths; or 2-vector [L1 L2], which is a shortcut for @(xy) L2-(L2-L1)*rssq(xy,2), i.e.
%       a linear mesh-size gradient from L2 at zenith to L1 at |xy| = 1 (usually horizon)
%
%     'projection','azimuthal' - See POLYPROJECTOR for options.
%
%     'offset',MESHSIZE/2,'cleantol',OFFSET/4 - Offset envelope by OFFSET, and simplify using
%       CLEANTOL. Ensures full interpolation coverage and uniform mesh close to edges.
%
%     'reduce', false - Include in V only points that would be used for interpolation with W,
%       i.e. return V(used,:) and W(:,used), where used = any(W,1). Also, if the generated mesh
%       has more vertices than the number of unique vectors S, the mesh will just be dumped and
%       replaced by unique(S), and the appropriate W.
%
%     'plot',false - meant for debugging
%
% See also: POLYPROJECTOR, INTERPMATRIX, GUESSOUTLINE, SOLARPOSITION.SUNPOSGRID
 
    MAX = 2.0;
    MIN = 0.8;

    if nargin == 0, test(); return; end

    % parse options / complete defaults
    [opt,varargin] = getflagoptions(varargin,{'-plot','-reduce'});

    opt.meshsize = getSimOption('angtol')*pi/180;
    opt.offset = opt.meshsize;
    opt.cleantol = opt.offset/4;
    opt.n = 6;
   
    [opt,varargin] = getpairedoptions(varargin,opt);
    parsestruct(opt,{'offset','n','cleantol','meshsize'},'-n','-r','-s','-p');
    
    % parse actual arguments
    isnice = @(x) isnumeric(x) && isreal(x) && all(isfinite(x(:)));
    switch numel(varargin)
        case 1
            P = varargin{1};
            if size(P,2)~=3, P = P'; end
            assert(ismatrix(P) && size(P,2) == 3,'Expecting 3-column S');
            P = P./rssq(P,2);
        case 2
            [az,el] = deal(varargin{:});
            az = az(:); el = el(:);
            assert(isnice(az) && isnice(el) && numel(az) == numel(el),'Bad azimuth-elevation');
            P = sph2cartV(az,el); %[Nt,3]
        otherwise
            error('Unexpected arguments')
    end
    assert(isnice(P),'S contains NaN/Inf/imaginary values');
    
    % Keep unique vectors only, idx works to return propper-sized W
    [P,~,idx] = unique(P,'rows');
    N = size(P,1);
        
    if N < 3
    % Can't generate a mesh with less than 3 points
        Q = P; 
        W = speye(N); W = W(idx,:);
        return; 
    end
    
    % PROVISIONAL?: Jigsaw crashes wiht a large mesh (1e5+)... and at some point it makes little 
    % sense to simplify a mesh rather than just starting with a coarse spherical tessellation
    N0 = 4*pi/(opt.meshsize^2*0.4);
    N0 = 2 + 10*4^ceil(log2((N0 - 2)/10)/2);
    if size(P,1) > N0
        Q = spherepoints(N0,'regular',true);
        if opt.plot, T = convhull(Q); end
    else    
        R = 4;
        ALPHA = R/2;

        P(N+1,:) = 0;  % add origin as last point

        shp = alphaShape(R*P,ALPHA*R);
        shp.HoleThreshold = shp.volume();
        T = shp.alphaTriangulation;
        T(~any(T > N,2),:) = [];        % remove tetrahedra without origin
        T = sort(T,2);                  % sort so that last vertex for each facet is origin

        T(:,4) = [];
        P(end,:) = [];

        TR = triangulation(T,P);
        E = TR.edges;
        toolong = dot(P(E(:,1),:),P(E(:,2),:),2) < cos(MAX*opt.meshsize);
        toobig = ismember(T(:,1:2),E(toolong,:),'rows') | ...
                 ismember(T(:,2:3),E(toolong,:),'rows') | ...
                 ismember(T(:,[1,3]),E(toolong,:),'rows');
        T(toobig,:) = [];

        warning_resetter = naptime('MATLAB:triangulation:PtsNotInTriWarnId'); %#ok<NASGU>

        [Q,T] = cleantrisurf(P,T,MIN*opt.meshsize,opt.meshsize,...
            'geom_feat',false,'optm_tria',false,'optm_div_',false,'mesh_rad2',10,'geom_seed',0);
     end

    [W,e] = interpmatrix(P,Q,'-sph','maxarc',opt.meshsize,'tol',opt.meshsize/16);
    W = blkdiag(W,speye(nnz(e)));
    Q = [Q;P(e,:)];

    if opt.reduce && size(Q,1) >= size(P,1)
    % More mesh vertices than original points: we'd be better off NOT interpolating
        Q = P; 
        W = speye(N); W = W(idx,:);
        return; 
    end
        
    if opt.reduce
    % Remove vertices that won't be used for interpolation
        used = any(W,1);
        Q = Q(used,:);
        W = W(:,used);

        if opt.plot
            used = full(sparse(find(used),1,1:nnz(used),numel(used),1));
            T = used(T);
            T(any(T==0,2),:) = [];
        end
    end
    W = W(idx,:); 

    if opt.plot        
        GUIfigure('sunposmesh'); clf(); hold on
        plothalfdome();
        trisurf(triangulation(T,Q),'facealpha',0.3,'edgecolor','b');
        scatter3(P(:,1),P(:,2),P(:,3),5,'r','fill');
        title(sprintf('%d points, reduced to %d vertices',size(P,1),size(Q,1)));
    end 
end

function test()
    rng(1);
    N = 1000;
    M = 10; 
    X = rand(N*M,3)*2-1 + (repmat(rand(M,3),N,1)*2-1)*5;
    X(:,3) = abs(X(:,3));
    X = X./rssq(X,2);
    tic()
    sunposmesh(X,'meshsize',4*pi/180,'-reduce','-plot');
    toc()
end