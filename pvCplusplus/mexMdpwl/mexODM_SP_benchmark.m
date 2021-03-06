
TOL = 1e-3;
setSimOption('RelTol',TOL);

PWLAPPROX = {'bisection','bijective','2nd',1/3}; % match to CPP for fair comparison
FULL = false;

% ODM = pickfile('*.samlib',Inf,'UI',1);
% ODM = cellfun(@pvlmod_SAMLibraryReader,ODM,'unif',0);
% ODM = cellfun(@checkODM,ODM,'unif',0);
% save('ODM_base_parameters.mat','ODM');

% Load ODM parameters & pick a module at random
load('ODM_base_parameters.mat','ODM');
% plot(cellfun(@(p) onediodepwlapprox(translateODM(p,800,40)),ODM));
% legend(arrayfun(@(j) sprintf('%d. %dc %s (%s)',j,ODM{j}.Ns,ODM{j}.material,ODM{j}.name),1:numel(ODM),'unif',0),'interpreter','none');
k = 1; % uncomment above to plot and make an informed choice
ODM = ODM{k};

% NOTE: unlike the (reduced & translated) parameter sets you had worked with before, these is
% now a complete parameter set at standard conditions, and has to be 'translated' to operating
% conditions by TRANSLATEODM.

% We want to simulate a shaded module, typically made up of 3 blocks of 20-24 cells each
Nb = ODM.Nbpd;                 % Number of cell-blocks
Ns = round(ODM.Ns/ODM.Nbpd);   % Number of cells per block

% Scale the parameters to reflect a single cell
CELL = scaleODM(ODM,1/ODM.Ns,1);  
CELL.Nbpd = 0;

% Generate IV-curve interpolation object
CIV = ModuleInterpolant(CELL,[],[],[],TOL,PWLAPPROX{:});

% Set absolute tolerances and limits based on the full-module's specs.
TOL = [ODM.Vmp_ref,ODM.Imp_ref,1]*TOL;
LIM = 1.5*[-ODM.Voc_ref ODM.Voc_ref -ODM.Isc_ref ODM.Isc_ref];

% Get the PWL of a bypass-diode, assume constant 40° heat-sink temperature for simplicity
D = bypassdiodefit('schottky');
v = bypassdiode2(D,2*ODM.Isc_ref,40,1e-6,'ambient');
[v,i] = pwlapprox(@(v) bypassdiode(D,v,40,1e-6,'ambient'),-ODM.Voc_ref,v,TOL,'bisection','euclidean','none');
diodeIV = mdpwl(-v,i);

%%
BLOCKFCN = {@(G,Tc) ODMseries(CELL,G,Tc,[],TOL,PWLAPPROX{:});
            @(G,Tc) ODMseries_interp(CIV,G,Tc,[],TOL);
            @(G,Tc) mexODMseries(CELL,G,Tc,[],TOL);
            @(G,Tc) mexVecODMseries(CELL,G,Tc,[],TOL)};
REPS = 100;
ModIV = {};
times = [];

for k = 1:REPS
    
    % Generate test irradiance & cell-temperature conditions for each cell (both Ns x Nb arrays)
    GBEAM = rand()*400+600;
    GSHADED = rand()*100+200;
    SHF = rand()*0.2;
    NOISE = 50; % as this value increases, performance advantage of mexVecODM... increases
    TA = 25;
    WS = 2.0;
    G = GSHADED + (rand([Ns,Nb]) > SHF)*GBEAM + randn([Ns,Nb])*NOISE;
    Tc = zeros(Ns,Nb);
    Tc(:) = pvl_sapmcelltemp(G(:),1000,-3.47,-0.0594,WS,TA,3);

    for i = 1:numel(BLOCKFCN)
        blckfcn = BLOCKFCN{i};

        blockIV = repmat(mdpwl(),Nb,1);
        tic()
        for j = 1:Nb
            % Get curveS for each (mismatched) string: onediodepwlapprox + mexaddseries
            blockIV(j) = blckfcn(G(:,j),Tc(:,j));
            if FULL, blockIV(j) = addparallel([blockIV(j),diodeIV],LIM); end
        end
        if FULL, ModIV{k,i} = addseries(blockIV); end
        times(k,i) = toc();
        % plot(ModIV)
    end
end
%%
GUIfigure('benchmark');
clf(); hold on;
boxplot(times/REPS);
set(gca,'yscale','log');
ylim([1e-5,1E-2]);
grid on;

% legend('onediodepwlapprox','ModuleInterpolant.getIVpp','mexVecODMpwlapproxAddseries');

function S = ODMseries(ODM,G,Tc,Lim,Tol,varargin)
% S = ODMSERIES(ODM,G,Tc,Lim,Tol) - Returns the single combined IV curve for a vector of identical
%   PV-elements (parameter structure ODM), operating at irradiances G and cell-temperatures Tc.
%   Loosely equivalent to the following calls:
%
%        P = translateODM(ODM,G,Tc);
%        C = arrayfun(@(p) onediodepwlapprox(p,Lim,Tol),P);
%        S = addseries(C,Lim,Tol);

    % narginchk(3,5);
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
        P = translateODM(ODM,G(j),Tc(j));                                       % translated parameters
        [XY{j}(:,1),XY{j}(:,2)] = onediodepwlapprox(P,Tol,Lim,varargin{:});    % individual IV curves
        XY{j}(:,1) = XY{j}(:,1)*W(j);                                           % .. w/ multiplicity
    end

    % Call MEX function, and wrap again into MDPWL object
    [xx,yy,sz] = mexaddseries(M,XY{:},Tol,Lim);
    S = mdpwl(xx(1:sz)', yy(1:sz)',0);
end

function S = ODMseries_interp(MIV,G,Tc,Lim,Tol)
% S = ODMSERIES_INTERP(M,G,Tc,Lim,Tol) - Same functionality as ODMSERIES, but using interpolation
%   from a grid of curves (MODULEINTERPOLANT object MIV).

    narginchk(3,5);
    ok = @(x) isvector(x) && all(isfinite(x)) && isreal(x);
    assert(ok(G) && ok(Tc) && numel(G) == numel(Tc),'Bad G,Tc');
    M = numel(G);
    if nargin < 4 || isempty(Lim), Lim = 1.5*[-MIV.Voc0*M MIV.Voc0*M -MIV.Isc0 MIV.Isc0]; end
    if nargin < 5 || isempty(Tol), Tol = [MIV.Vmp0,MIV.Imp0,1]*getSimoption('RelTol'); end
    
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
    
    % Get IV curves from interpolant, apply multiplicity
    P = MIV.getIVpp(G,Tc);
    XY = arrayfun(@(p,w) [p.x*w,p.y],P,W,'unif',0);

    % Call MEX function, and wrap again into MDPWL object
    [xx,yy,sz] = mexaddseries(M,XY{:},Tol,Lim);
    S = mdpwl(xx(1:sz)', yy(1:sz)',0);
end

function S = mexVecODMseries(ODM,G,Tc,Lim,Tol)
% S = MEXODMSERIES(ODM,G,Tc,Lim,Tol) - Same functionality as ODMSERIES, 
%   but using mexVecODMpwlapproxAddseries

    narginchk(3,5);
    ok = @(x) isvector(x) && all(isfinite(x)) && isreal(x);
    assert(ok(G) && ok(Tc) && numel(G) == numel(Tc),'Bad G,Tc');
    M = numel(G);
    if nargin < 4 || isempty(Lim), Lim = 1.5*[-ODM.Voc_ref*M ODM.Voc_ref*M -ODM.Isc_ref ODM.Isc_ref]; end
    if nargin < 5 || isempty(Tol), Tol = [ODM.Vmp_ref,ODM.Imp_ref,1]*getSimoption('RelTol'); end
    
    % In practice, many modules/cells will have very similar irradiance, so it makes sense to
    % identify them and treat them as a single (scaled) IV curve
    GVARTOL = sqrt(Tol(end));  % a function of the slope at MPP, but this should do for now
    [~,ic,ia] = uniquetol(G,GVARTOL);
    if numel(ic) < M
        W = accumarray(ia,1);
        G = accumarray(ia,G,[],@mean);
        Tc = accumarray(ia,Tc,[],@mean);
    else
        W = ones(M,1);
    end
    
    % combination of vector handling and addseries
    P = translateODM(ODM,G,Tc);

    [xx, yy, sz] = mexVecODMpwlapproxAddseries(P(:)', Tol(:)', Lim(:)', W(:));
    S = mdpwl(xx(1:sz)', yy(1:sz)',0);
end

function S = mexODMseries(ODM,G,Tc,Lim,Tol)
% S = MEXODMSERIES(ODM,G,Tc,Lim,Tol) - Same functionality as ODMSERIES, 
%   but using mexVecODMpwlapproxAddseries

    narginchk(3,5);
    ok = @(x) isvector(x) && all(isfinite(x)) && isreal(x);
    assert(ok(G) && ok(Tc) && numel(G) == numel(Tc),'Bad G,Tc');
    M = numel(G);
    if nargin < 4 || isempty(Lim), Lim = 1.5*[-ODM.Voc_ref*M ODM.Voc_ref*M -ODM.Isc_ref ODM.Isc_ref]; end
    if nargin < 5 || isempty(Tol), Tol = [ODM.Vmp_ref,ODM.Imp_ref,1]*getSimoption('RelTol'); end
    
    % In practice, many modules/cells will have very similar irradiance, so it makes sense to
    % identify them and treat them as a single (scaled) IV curve
    GVARTOL = sqrt(Tol(end));  % a function of the slope at MPP, but this should do for now
    [~,ic,ia] = uniquetol(G,GVARTOL);
    if numel(ic) < M
        W = accumarray(ia,1);
        G = accumarray(ia,G,[],@mean);
        Tc = accumarray(ia,Tc,[],@mean);
        M = numel(ic);
    else
        W = ones(M,1);
    end
    
    XY = cell(1,M);
    for j = 1:M
        P = translateODM(ODM,G(j),Tc(j));                        % translated parameters
        [XY{j}(:,1),XY{j}(:,2)] = mexODMpwlapprox(P,Tol,Lim);    % individual IV curves
        XY{j}(:,1) = XY{j}(:,1)*W(j);                            % .. w/ multiplicity
    end

    % Call MEX function, and wrap again into MDPWL object
    [xx,yy,sz] = mexaddseries(M,XY{:},Tol,Lim);
    S = mdpwl(xx(1:sz)', yy(1:sz)',0);
end
