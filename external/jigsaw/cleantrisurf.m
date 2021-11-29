function [V,T] = cleantrisurf(V,T,varargin)
% [V,T] = CLEANTRISURF(V,T,HMIN,HMAX,...) - JIGSAW wrapper for surface mesh simplification.
%   See JIGSAW (scarce) documentation for options.
%
% NOTES: 

    if nargin == 0, testjigsaw(0); return; end
    
    % HFUN_PTS = 32;
    
    opts.hfun_hmin = [];
    opts.hfun_hmax = [];

    opts.init_near = []; % 1e-8;
    opts.hfun_scal = 'absolute'; % 'relative'
    opts.mesh_kern = 'delaunay'; % 'delfront' (hangs with big hmax)
    
    opts.geom_seed = []; % 8
    
    opts.geom_phi1 = []; % 6 (undocumented)
    opts.geom_phi2 = []; % 6 (undocumented)
    opts.geom_eta1 = []; % 45 angle that makes a "sharp" feature
    opts.geom_eta2 = []; % 45 
    opts.geom_feat = []; % false;
    
    opts.mesh_top1 = true;
    opts.mesh_top2 = true;
    opts.mesh_eps1 = []; % 0.33;
    opts.mesh_eps2 = []; % 0.33;
    opts.mesh_rad2 = 1.5; % 1.05;
    opts.mesh_off2 = []; % 0.9;
    opts.mesh_snk2 = []; % 0.2;
    
    % opts.mesh_rad3 = 2.05;
    % opts.mesh_off3 = 1.1;
    % opts.mesh_snk3 = 0.33;
    % opts.mesh_vol3 = 0;
    
    opts.optm_kern = 'odt+dqdx'; % 'cvt+dqdx'
    opts.optm_tria = true;
    opts.optm_dual = false;
    opts.optm_zip_ = true;
    opts.optm_div_ = true;
    opts.optm_iter = 16;
    opts.optm_qtol = 1e-4;
    
    opts.verbosity = -1;
    % opts.hfun = [];
    
    opts = getpairedoptions(varargin,opts,'dealrest',2);
    % hfun = opts.hfun; opts.hfun = [];
    fld = fieldnames(opts);
    opts = rmfield(opts,fld(cellfun(@isempty,struct2cell(opts))));
    
    opts.mesh_dims = 2 ;
    
    validateattributes(V,{'numeric'},{'real','finite','size',[NaN,3]},'','V');
    validateattributes(T,{'numeric'},{'real','integer','positive','<=',size(V,1),'size',[NaN,3]},'','T');
    % parsestruct(opts,{'hfun_hmin','hfun_hmax'},'numeric','real','scalar','positive');
    % validateattributes(tol,{'numeric'},{'real','scalar','positive','<',1},'','tol');
    
    % if ~isempty(hfun)
    %     validateattributes(hfun,{'numeric','function_handle','griddedInterpolant'},...
    %         {'scalar'},'','hfun');
    % else
    %     hfun = []; 
    % end
    
    initjigsaw();
    % global JIGSAW_EDGE2_TAG
    
    rootpath = fileparts(which('jigsaw')) ;
    
    TR = triangulation(T,V);
    C = TR.freeBoundary;

    geom.point.coord = [V,zeros(size(V,1),1)];
    geom.tria3.index = [T,zeros(size(T,1),1)];
    geom.edge2.index = [C,-ones(size(C,1),1)];
    geom.mshID = 'EUCLIDEAN-MESH';
    % geom.mshID = 'ELLIPSOID-MESH';
    % geom.radii = [1,1,1];

    name = 'cleantrisurf';
    opts.geom_file = fullfile(rootpath,'files',[name,'.msh']) ; % domain file
    opts.jcfg_file = fullfile(rootpath,'cache',[name,'.jig']) ; % config file
    opts.mesh_file = fullfile(rootpath,'cache',[name,'.msh']) ; % output file
    
    savemsh(opts.geom_file,geom);
    
    % FIX: jigsaw is not working with hfun_files, no idea why.
    
    % if ~isempty(hfun) && isnumeric(hfun)
    %     H = hfun;
    %     hfun = @(x,y) repmat(H,size(x));
    %     HFUN_PTS = 2;
    % end
    % if ~isempty(hfun)
    %     try
    %         if isa(hfun,'griddedInterpolant')
    %             hmat.point.coord = hfun.GridVectors;
    %             hmat.value = hfun.Values;
    %         else
    %             [lo,hi] = bounds(V,1);
    %             hmat.point.coord{1} = linspace(lo(1),hi(1),HFUN_PTS);
    %             hmat.point.coord{2} = linspace(lo(2),hi(2),HFUN_PTS);
    %             [X,Y] = meshgrid(hmat.point.coord{1:2});
    %             hmat.value = hfun(X,Y);
    %             validateattributes(hmat.value,{'numeric'},{'real','size',size(X)},'','hfun(X,Y)');
    %         end
    %     catch
    %         error('Invalid distance function handle');
    %     end        
    %     hmat.mshID = 'EUCLIDEAN-GRID';
    % 
    %     opts.hfun_file = fullfile(rootpath,'cache',[name,'-hfun.msh']) ;
    %     savemsh(opts.hfun_file,hmat) ;
    % end

    if opts.verbosity >= 0
        mesh = jigsaw(opts);
    else
        [~,mesh] = evalc('jigsaw(opts)');
    end
    assert(isfield(mesh,'point'),'jigsaw failed, check your options'); 
    
    V = mesh.point.coord(:,1:3);
    T = mesh.tria3.index(:,1:3);
    
    TR = triangulation(T,V);
    n = TR.faceNormal((1:size(T,1))');
    facedown = n(:,3) < 0;
    T(facedown,:) = fliplr(T(facedown,:));
        
    endjigsaw();
end

