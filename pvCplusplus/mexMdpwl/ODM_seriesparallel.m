
TOL = 1e-4;
setSimOption('RelTol',TOL);

% Load ODM parameters & pick a module at random
addpath(relativepath('../../../tests/ODM',fileparts(mfilename('fullpath'))));
ODM = ODM_testbase();
k = randi(numel(ODM),1); 
ODM = ODM{k};

% NOTE: unlike the (reduced & translated) parameter sets you had worked with before, these is
% now a complete parameter set at standard conditions, and has to be 'translated' to operating
% conditions by TRANSLATEODM.

% We want to simulate a shaded module, typically made up of 3 blocks of 20-24 cells each
Nb = ODM.Nbpd;          % Number of cell-blocks
Ns = round(ODM.Ns/ODM.Nbpd);   % Number of cells per block

% Scale the parameters to reflect a single cell
CELL = scaleODM(ODM,1/ODM.Ns,1);  

% Set absolute tolerances and limits based on the full-module's specs.
TOL = [ODM.Vmp_ref,ODM.Imp_ref,1]*TOL;
LIM = 1.5*[-ODM.Voc_ref ODM.Voc_ref -ODM.Isc_ref ODM.Isc_ref];

% Get the PWL of a bypass-diode, assume constant 40° heat-sink temperature for simplicity
D = bypassdiodefit('schottky');
v = bypassdiode2(D,2*ODM.Isc_ref,40,1e-6,'ambient');
[v,i] = pwlapprox(@(v) bypassdiode(D,v,40,1e-6,'ambient'),-ODM.Voc_ref,v,TOL);
diodeIV = mdpwl(-v,i);

% Generate test irradiance & cell-temperature conditions for each cell (both Ns x Nb arrays)
TA = 25;
WS = 2.0;
G = 200 + (rand([Ns,Nb]) > 0.05)*800 + randn([Ns,Nb])*2.0;
Tc = zeros(Ns,Nb);
Tc(:) = pvl_sapmcelltemp(G(:),1000,-3.47,-0.0594,WS,TA,3);


blockIV = repmat(mdpwl(),Nb,1);
tic()
for j = 1:Nb
    % Get the curve for each (mismatched) string at once
    blockIV(j) = future_mexODMseries(CELL,G(:,j),Tc(:,j),[],TOL);
end
toc()
clf()
plot(blockIV);

% The addition of each cell-block with its bypass diode; then in series with other blocks/
% modules to form strings; and in parallel with other strings to form an array-IV-curve would
% be the second step.
tic()
for j = 1:Nb
    blockIV(j) = addparallel([blockIV(j),diodeIV],LIM);
end
ModIV = addseries(blockIV);
toc()
plot(blockIV);
plot(ModIV)


function S = future_mexODMseries(ODM,G,Tc,Lim,Tol)
% S = ODMSERIES(ODM,G,Tc,Lim,Tol) - Returns the single combined IV curve for a vector of identical
%   PV-elements (parameter structure ODM), operating at irradiances G and cell-temperatures Tc.
%   Loosely equivalent to the following calls:
%
%        P = translateODM(ODM,G,Tc);
%        C = arrayfun(@(p) onediodepwlapprox(p,Lim,Tol),P);
%        S = addseries(C,Lim,Tol);
%
% This would be the next function to port into C++. Note that it includes a call to ODMTRANSLATE.
% We could go around it if I pass to you an array of structures, or even a matrix, with the
% parameters already 'translated'

    narginchk(3,5);
    ok = @(x) isvector(x) && all(isfinite(x)) && isreal(x);
    assert(ok(G) && ok(Tc) && numel(G) == numel(Tc),'Bad G,Tc');
    M = numel(G);
    if nargin < 4 || isempty(Lim), Lim = 1.5*[-ODM.Voc_ref*M ODM.Voc_ref*M -ODM.Isc_ref ODM.Isc_ref]; end
    if nargin < 5 || isempty(Tol), Tol = [ODM.Vmp_ref,ODM.Imp_ref,1]*getSimoption('RelTol'); end
    
    % In practice, many modules/cells will have very similar irradiance, so it makes sense to
    % identify them and treat them as a single IV curve.
    GVARTOL = sqrt(Tol(end));  % a function of the slope at MPP, but this should do for now
    [~,ic,ia] = uniquetol(G,GVARTOL);
    if numel(ic) < M
        W = accumarray(ia,1);
        G = accumarray(ia,G,[],@mean);
        Tc = accumarray(ia,Tc,[],@mean);
        M = numel(ic);
    else
        % M = M;
        W = ones(M,1);
    end
    XY = cell(1,M);
    for j = 1:M
        P = translateODM(ODM,G(j),Tc(j));                           % translated parameters
        [XY{j}(:,1),XY{j}(:,2)] = onediodepwlapprox(P,Tol,Lim);     % individual IV curves
        XY{j}(:,1) = XY{j}(:,1)*W(j);                               % .. w/ multiplicity
    end

    % Call MEX function, and wrap again into MDPWL object
    [xx,yy,sz] = mexaddseries(M,XY{:},Tol,Lim);
    S = mdpwl(xx(1:sz)', yy(1:sz)',0);
end