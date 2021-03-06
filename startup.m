function startup(try_to_compile)
if nargin < 1, try_to_compile = false; end

%% Add paths

BASEDIR = fileparts(mfilename('fullpath'));
HERE = pwd();
lastwill = onCleanup(@() cd(HERE));

cd(BASEDIR);
addpath(genpath('.'));
rmpath(genpath('./.git'));
rmpath(genpath('./tests/ODM'));

rmpath(genpath('./pvCplusplus'));
addpath('./pvCplusplus');
addpath('./pvCplusplus/mexMdpwl');

rmpath('./matlab-pvmeteo/PVLIB/Example Data');
rmpath('./external/clipper/cpp');
rmpath(genpath('./external/jigsaw/jigsaw-matlab'));

%% Compile CPP libraries

iscompiled = @(x) endsWith(which(x),{'mexa64','mexw64','mexa32','mexw32'});

% Compile C++ One-Diode-Model and Piece-Wise-Linear electrical solver

if ~all(cellfun(iscompiled,{'mexaddseries','mexaddparallel','mexVecODMpwlapproxAddseries','mexODMpwlapprox'}))
    if try_to_compile
        try
            cd(fullfile(BASEDIR, 'pvCplusplus/mexMdpwl'));
            commonBuildScript()
            addSeriesBuildScript()
            addParallelBuildScript()
            ODMpwlapproxBuildScript()
        catch
            warning('Failed to compile C++ ODM and PWL solver')
        end
    else
        warning('C++ ODM and PWL solver still need to be compiled')
    end
end

% Compile CLIPPER polygon-clipping library/mex wrapper (see `/external/clipper/matlab/README.txt`)
if ~iscompiled('clipper')
    if try_to_compile
        try
            cd(fullfile(BASEDIR,'external/clipper/matlab'));
            cd('./external/clipper/matlab');
            mex '-D__int64=__int64_t' '-I../cpp' '../cpp/clipper.cpp' 'mexclipper.cpp'
        catch
            warning('Failed to compile clipper')
        end
    else
        warning('Clipper still needs to be compiled')
    end
end

% Compile JIGSAW mesh refining utility (see `/external/jigsaw/jigsaw-matlab/README.md`)

cd(fullfile(BASEDIR,'external/jigsaw'));
initjigsaw();
try
    evalc('example(9,0)');
catch
    if try_to_compile
        try
            compile();
        catch
            warning('Failed to compile jigsaw')
        end
    else
        warning('Jigsaw still needs to be compiled')
    end
end
endjigsaw();

end
