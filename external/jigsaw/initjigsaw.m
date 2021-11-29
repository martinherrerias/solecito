function initjigsaw()
% just a gateway for initjig

    filename = mfilename('fullpath') ;
    filepath = fileparts( filename ) ;
    
    addpath([filepath, '/jigsaw-matlab']) ;
    initjig();
end



