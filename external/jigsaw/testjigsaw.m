function testjigsaw(verbosity)
% internal JIGSAW test for surface mesh simplification

    if nargin < 1, verbosity = -1; end

    jexename = fileparts(mfilename('fullpath'));
    if ispc()
        jexename = [jexename,'\jigsaw-matlab\external\jigsaw\bin\jigsaw.exe'];
    else
        jexename = [jexename,'/jigsaw-matlab/external/jigsaw/bin/jigsaw'];
    end
    assert(isfile(jexename), ['JIGSAW''s executable not found -- ', ...
        'See jigsaw-matlab/compile.m for additional detail.']);
    
    blob = polygon(sort(rand(10,1))*360,movmean(rand(10,1),2),'pol');
    blob.x = blob.x*100+100;
    blob.y = blob.y*100+100;
    blob = polygon(interpn(blob.x([1:end,1]),3,'makima'),interpn(blob.y([1:end,1]),3,'makima'));
    % c = centroid(blob);
    
    [x,y] = deal(rand(10000,1)*200,rand(10000,1)*200);
    ix = inpolygon(x,y,blob.x,blob.y);
    x = x(ix); y = y(ix);
    
    z = cos(hypot(2*x,y)/20)*5;
    C = boundary(x,y,0.2);
    C(:,2) = circshift(C,-1);
    T = delaunayTriangulation([x,y],C(1:end-1,:));
    T = T.ConnectivityList(T.isInterior(),:);
    
    % hfun = @(x,y) 1 - 0.8*exp(-((x-c(1)).^2 + (y-c(2)).^2)/5000) ;
    
    [V,T] = cleantrisurf([x,y,z],T,...
        'hfun_hmin',0.1,'hfun_hmax',10,'verbosity',verbosity); % 'mesh_eps1',0.02,'mesh_eps2',0.02,
    
    if verbosity >= 0
        GUIfigure('jigsaw','Jigsaw test','2:1'); clf(); hold on;
        trisurf(T,V(:,1),V(:,2),V(:,3));
        plot3(x,y,z,'r.')
        view(45,45)
    else
       % [log,mesh] = evalc('jigsaw(opts)');
       % assert(~isempty(mesh),log);
    end
end
