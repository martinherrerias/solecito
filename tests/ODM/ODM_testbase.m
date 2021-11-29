function [ODMs,MODULES] = ODM_testbase(tol)

warning_resetter = naptime({'checkODM:assumption','checkODM:nomchk'}); %#ok<NASGU>
if nargin < 1 || isempty(tol),tol = getSimOption('RelTol'); end

wd = pwd(); workspace_resetter = onCleanup(@() cd(wd));
cd(fileparts(mfilename('fullpath')));
files = pickfile('*.samlib',Inf);

setSimOption('RelTol',tol);
setSimOption('MaxIter',200);

ODMs = cell(numel(files),1);
for j = numel(ODMs):-1:1
    ODMs{j} = pvlmod_SAMLibraryReader(files{j});
    ODMs{j} = checkODM(ODMs{j});
end
if nargout < 2, return; end

subdir = sprintf('%0.0e',tol);
if ~isfolder(subdir), mkdir(subdir); end
cd(subdir);
for j = numel(ODMs):-1:1
    if isempty(dir([ODMs{j}.name '.modint']))
    % Generate MODULEINTERPOLANT objects
       ModIVint = ModuleInterpolant(ODMs{j});
       save([ODMs{j}.name '.modint'],'ModIVint');
    else
    % ... or load precalculated
       load([ODMs{j}.name '.modint'],'-mat','ModIVint');
    end
    %pp = ModIVint.getIVpp(1000,40);
    %pp.clip([0 Inf 0 Inf]);
    %plot(pp);
    MODULES(j) = ModIVint;
end
