function test()
% internal JIGSAW test for surface mesh simplification

    initjig ;
    rootpath = fileparts(which('initjig')) ;
    name = 'test';


    [x,y] = deal(rand(10000,1)*200,rand(10000,1)*200);
    z = 4*cos(sqrt((3*x/60-2).^2 + (2*y/60+1).^2));
    C = boundary(x,y,0.2);
    C(:,2) = circshift(C,-1);
    T = delaunayTriangulation([x,y],C(1:end-1,:));
    T = T.ConnectivityList(T.isInterior(),:);
    
    geom.point.coord = [x,y,z,zeros(numel(x),1)];
    geom.tria3.index = [T,zeros(size(T,1),1)];
    geom.mshID = 'EUCLIDEAN-MESH';

    opts.geom_file = fullfile(rootpath,'files',[name,'.msh']) ;  % domain file
    opts.jcfg_file = fullfile(rootpath,'cache',[name,'.jig']) ; % config file
    opts.mesh_file = fullfile(rootpath,'cache',[name,'.msh']) ; % output file

    opts.mesh_dims = 2 ;
    opts.hfun_hmax = 10.0 ;
    opts.hfun_hmin = 0.5 ;
    opts.hfun_scal = 'absolute';
    opts.mesh_kern = 'delaunay'; % / 'delfront';
    opts.mesh_top1 = true ;
    opts.geom_feat = true ;
    % opts.mesh_eps1 = 0.1/opts.hfun_hmax;
    opts.mesh_eps2 = 0.2/opts.hfun_hmin;

    savemsh(opts.geom_file,geom);

    %    [XPOS,YPOS] = meshgrid(xpos,ypos) ;
    %     hfun =-.4*exp(-.1*(XPOS-4.5).^2 ...
    %                   -.1*(YPOS-4.5).^2 ...
    %             ) + .6 ;
    % 
    %     hmat.mshID = 'EUCLIDEAN-GRID';
    %     hmat.point.coord{1} = xpos ;
    %     hmat.point.coord{2} = ypos ;
    %     hmat.value = hfun ;
    %     savemsh(opts.hfun_file,hmat) ;


    GUIfigure('jigsaw','Jigsaw test','2:1'); clf(); hold on;
    % drawmesh(geom);

    [~,mesh] = evalc('jigsaw(opts)');
    drawmesh(mesh);
    plot3(x,y,z,'r.')
    view(45,45)
end
