function endjigsaw()
% revert INITJIG, to avoid possible conficts with other libraries

    filename = mfilename('fullpath') ;
    filepath = fileparts( filename ) ;
    
    rmpath( [filepath, '/jigsaw-matlab'] ) ;
    rmpath( genpath( [filepath, '/jigsaw-matlab/main'] )) ;
    rmpath( genpath( [filepath, '/jigsaw-matlab/tools'] )) ;
    rmpath( genpath( [filepath, '/jigsaw-matlab/parse'] )) ;

    clear GLOBAL -REGEXP JIGSAW.*
end



