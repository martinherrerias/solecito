function [POA,EB]= effectiveirradiance(POA,MD,Trck,pvModType)
% [POA*,EB] = EFFECTIVEIRRADIANCE(POA,MD,SP,TRCK) - Interpolate the results of POAIRRADIANCE to 
%   return module-average diffuse irradiance values, and apply modifying factors (spectral 
%   correction, soiling, degradation...) to get diffuse- and beam-effective irradiances.
%   Keep track of the chain of losses on a (mount averaged) Energy-Balance (EB) structure, which
%   also includes approximations for the electrical shading losses, based on linear shading
%   factors and number of shaded cell-blocks.
%
% [~,EB] = EFFECTIVEIRRADIANCE(POA*,MD,SP,TRCK) - recalculate the Energy-Balance (EB) structure
%   from a previous run (already module-averaged and corrected POA*).
%
% Input: 
% MD - structure with Nt-vector fields DHI (Diffuse Horizontal-), and GHI (Global Horizontal-)
%       Irradiance [W/m²]; and optionally AMa (pressure-corrected air mass) and tpw (precipitable 
%       water content) [kg/m²] for spectral correction. MD.soiling and MD.degradation are optional
%       scalar of Nt-vector fields that will be applied as linear irradiance loss factors 
%       (1 = no irradiance, default = 0). For cell/module specific degradation consider using
%       POA.Defects inside PVARRAYSOLVER.
%
% TRCK.geom - pvArea object of depth 4 representing the tracker and its sub-divisions in modules,
%           cell-blocks, and cells.
% TRCK.analysedpoints - [2,Npt] coordinates of shading-analysis nodes in each mount.
%
% pvModType - module material, used by pvl_FSspeccorr
%
% OUTPUT:
% 	POA*.{ Bpoa0, Dpoa, Bpoa_IAM, Dpoa_IAM, Nsb } fields untouched from POA (see POAIRRADIANCE)
% 	POA*.Bpoa - Nt·Nu* array with intensity of direct irradiance after diffuse-shading, soiling, 
%       spectral- and incidence angle correction. SHADING NOT INCLUDED
%   POA*.ShF - Nt·Nu·Nm array with fractions of shade-coverage for each module (1..Nm) in every
%       tracker (1:Nu), at all timesteps. ("Linear" shading, 1 = full shade, 0 = no shade).
% 	POA*.Dpoa - Nt·Nu·Nm array with intensity of diffuse irradiance components after shading, 
%       soiling, spectral- and incidence angle correction. Averaged for every MODULE in the array.
% 	POA*.Dstd - Nt·Nu·Nm array with intra-module diffuse irradiance variability (std. dev.) to
%       be used for small mismatch correction (e.g. Bucciarelli, 1979).
%   The Energy-Balance (EB) structure contains Nt-vector fields for complete-plant (all modules
%   in all analysed mounts) averaged values of irradiance and irradiance losses.
%
%   EB.area - (scalar) m² of modules, assuming a fully-connected array
% 	EB.GHI - Global horizontal irradiance over active module area [W/m²]
%   EB.bias - Corrected bias due to far-horizon & projection (BCF in POAIRRADIANCE)
% 	EB.Bpoa0 - POA beam irradiance over active module area, without shading [W/m²]
% 	EB.Dpoa0 - POA diffuse irradiance over active module area, without shading [W/m²]
%   EB.IAMb - Absortivity Losses due to Incidence-Angle of beam irradiance [W/m²]
%   EB.IAMd - Absortivity Losses due to Incidence-Angle of diffuse irradiance [W/m²]
%   EB.ShLb - Linear beam irradiance losses (no mismatch effects)
% 	EB.ShLd - Linear diffuse irradiance losses (no mismatch effects)
%   EB.ShLaE - Approx. Electrical shading losses [Martínez-Moreno et al. (2010)], central inverter
%   EB.ShLaU - Approx. Electrical shading losses, using mount == string inverters
%   EB.ShLaW - Worst case Electrical shading losses
%   EB.Gpoa - Global POA irradiance using linear shading factors [W/m²]	
% 	EB.spectral - Spectral Correction Loss [W/m²]
% 	EB.soiling - Soiling Loss [W]
% 	EB.De - Effective POA Diffuse Irradiance, after IAM, spectral & soiling losses [W]
%   EB.Be - Effective POA Beam Irradiance, after IAM, spectral & soiling losses NOT SHADED [W]
%   EB.Ge - Effective POA Global Irradiance, after spectral & soiling losses, linear shading [W]    
%
% Martínez-Moreno, F., Muñoz, J., Lorenzo, E., 2010. "Experimental model to estimate shading 
% losses on PV arrays". Solar Energy Materials and Solar Cells 94, 2298–2303.
% L. L. Bucciarelli, Power loss in photovoltaic arrays due to mismatch in cell characteristics
%   Solar Energy, vol. 23, no. 4, pp. 277–288, Jan. 1979, doi: 10.1016/0038-092X(79)90121-X.
%
% See also: POAIRRADIANCE, CHECKIAM, GUIIRRTRANS

    narginchk(3,4);
    if nargin < 4, pvModType = ''; end
    
    parsestruct(Trck,{'analysedpoints'},'-n','-r','size',[2,NaN])
    parsestruct(Trck,{'geom'},'class','pvArea',@(x) x.depth == 4);    
    parsestruct(MD,{'GHI','DHI'},'opt',{'AMa','tpw'},'-n','-r','-v','-e');
    
    [Nm,Nb,~] = Trck.geom.getdims();    % module positions, cell-blocks per module 
    modulearea = Trck.geom.elements(1).area;
    
    Nt = size(MD.GHI,1); % time-steps
    
    isshaded = all(isfield(POA,{'ShF','Dpoa','Nsb'}));
    iseffective = all(isfield(POA,{'Dstd','Bpoa'}));
        
    if ~iseffective && isshaded
    % Interpolate from analysis-points to modules
    
        [POA.Dpoa,POA.Dstd] = nodalpoints2modules(POA.Dpoa,Trck.geom,Trck.analysedpoints);
        if size(POA.CS,3) ~= size(POA.Dpoa,3)
        % FIX! it could happen that Nm == Np, but still CS was view-point based, then we are
        %   in trouble.
            POA.CS = nodalpoints2modules(POA.CS,Trck.geom,Trck.analysedpoints);
        end
    end
    
    % At this point all fields should be singleton-expansion compatible [Nt,Nu*,Nm*]
    parsestruct(POA,{'Bpoa0','Dpoa0','Bpoa_IAM','Dpoa_IAM','bias','CS'},...
        'opt',{'Bpoa','Dstd','ShF','Nsb','Dpoa'},'-n','-r','-c','size',[Nt,NaN,NaN]);
    
    if isshaded
        Nu = size(POA.Dpoa,2);  % number of analysed trackers
        assert(size(POA.ShF,3) == Nm,'Inconsistent module numbers');
    else
        Nu = size(POA.Dpoa_IAM,2);
    end
    
    % Define integrating functions for Energy-Balance reporting:
    EB.area = Nu*Nm*modulearea;
    mountmean = @(I) mean(I,2);     % assume all modules are connected
    modmean = @(I) mean(I,2:3);
        
    % Energy Balance: start from horizontal irradiance...
	EB.GHI = mountmean(repmat(MD.GHI,1,Nu));
    
    % Record bias due to far-horizon & projection
    EB.bias = (1-POA.bias).*MD.DHI;
    EB.bias(MD.GHI == 0) = 0;

    % Unshaded components (IAM = 1)
	EB.Bpoa0 = mountmean(POA.Bpoa0); % [Nt,Nu]
	EB.Dpoa0 = mountmean(POA.Dpoa0);
        
    % Approx IAM losses (before shading). 
    EB.IAMb = mountmean(POA.Bpoa0 - POA.Bpoa_IAM);
    EB.IAMd = mountmean(POA.Dpoa0 - POA.Dpoa_IAM);
    
    Gpoa0 = mean(POA.Bpoa_IAM + POA.Dpoa_IAM,2);
    dark = Gpoa0 <= 0;
    Gpoa0(dark) = 0;
    
    % Soiling
    if isfield(MD,'soiling')
        fsoiling = MD.soiling;
    else
        fsoiling = 0;
        % warning('Applying no soiling correction');
    end
	EB.soiling = modmean(fsoiling.*Gpoa0);

    % Spectral correction
    if ~all(isfield(MD,{'tpw','AMa'}))
        pvModType = '';
        warning('AMa and tpw not available for spectral correction');
    end
    pvModType = checkmodType(pvModType); % includes warning if isempty(pvModType)
    if ~isempty(pvModType)
        fspectral = 1 - pvl_FSspeccorr(MD.tpw/10,MD.AMa,pvModType);
        fspectral(~isfinite(MD.AMa)) = 0;
    else
        fspectral = zeros(Nt,1);
    end
    EB.spectral = modmean(fspectral.*Gpoa0);
    
    % Intra-module mismatch due to diffuse (modified Bucciarelli, 1979).
    % TODO: currently only informative (printed to EB), should be carried over and resolved by
    % PVARRAYSOLVER. Problem: even though mismatch at PMP scales with (std(G)/G)² the IV curve
    % for V < Vmp can be affected by O(std(G)/G). See T20210228_Mismatch_again.m on log.
    Ns = prod(Trck.geom.dims(2:3));
    c = 18.7 - 0.0054*Gpoa0; % empirical fit (mSi,pSi), ballpark for any module technology
    fmismatch = (1-1/Ns).*(c/2+1).*(POA.Dstd./Gpoa0).^2;
    fmismatch = min(POA.Dpoa,fmismatch);
    EB.bucciarelli = modmean(fmismatch);
    
    % Linear diffuse and beam irradiance losses
    if ~iseffective
        EB.ShLd = modmean(POA.Dpoa_IAM - POA.Dpoa);
    else
        EB.ShLd = modmean(POA.Dpoa_IAM - POA.Dpoa./(1-fspectral-fsoiling));
    end
    EB.ShLb = modmean(POA.ShF.*POA.Bpoa_IAM);
    EB.ShLd(dark) = 0;
    EB.ShLb(dark) = 0;
    
    if ~isshaded, POA.Dpoa = POA.Dpoa_IAM; end
    if ~iseffective
        POA.Dpoa = POA.Dpoa.*(1-fspectral-fsoiling);
        POA.Bpoa = POA.Bpoa_IAM.*(1-fspectral-fsoiling);
    end
        
    % Approx. Electrical shading losses (ShLa..) according to Martínez-Moreno et al. (2010)
    BpCS = mean((1-POA.ShF).*POA.Bpoa + POA.CS,3);    % (1-Fgs)(B + CS) [Nt,Nu]
    DmCS = mean(POA.Dpoa - POA.CS,3);                 % SKY + GND       [Nt,Nu]

        % U. Unit mounts form individual arrays, e.g. mount/string inverters
        fe = 1-sum(POA.Nsb,3)/(Nm*Nb+1);  % trck missmatch factor (1-Fes)/(1-Fgs)
        EB.ShLaU = Gpoa0-mean(BpCS.*fe + DmCS,2);
        EB.ShLaU(dark) = 0;

        % E. Everything is connected into a single array (central inverter)
        fe = 1-sum(POA.Nsb,2:3)/(Nu*Nm*Nb+1);  % whole array missmatch factor (1-Fes)/(1-Fgs)
        EB.ShLaE = Gpoa0-mean(BpCS.*fe + DmCS,2);
        EB.ShLaE(dark) = 0;
    
    % Worst Case
    EB.ShLaW = Gpoa0 - min(POA.Dpoa + (POA.Nsb == 0).*POA.Bpoa,[],2:3);

    % Module-level Gpoa using shaded diffuse irradiance and linear-shaded beam ($)
    Gpoa = (1-POA.ShF).*POA.Bpoa + POA.Dpoa;
    EB.Gpoa = modmean(Gpoa);

    % Global equivalent components
	EB.De = modmean(POA.Dpoa);
    EB.Be = mountmean(POA.Bpoa);
    EB.Ge = modmean(POA.Dpoa + (1-POA.ShF).*POA.Bpoa);
end

function [M,S] = nodalpoints2modules(Q,geom,pts)
% Interpolate variable Q [Nt,Ntr,Npt] evaluated at nodal-points PTS, 
% at the center of every module in GEOM, returning an [Nt,Ntr,Nm]
% interpolated array M.

    % Get cell centers
    sz = geom.dims;
    C = zeros([sz,2]);
    for i = 1:sz(1)
        for j = 1:sz(2)
            for k = 1:sz(3)
                C(i,j,k,1) = mean(geom.elements(i).elements(j).elements(k).border.x);
                C(i,j,k,2) = mean(geom.elements(i).elements(j).elements(k).border.y);
            end
        end
    end
    
    Nm = size(C,1);
    Ns = sz(2)*sz(3);
    [Nt,Ntr,Npt] = size(Q);
    assert(Npt == size(pts,2),'Dpoa size does not match node-point coordinates');
    
    C = reshape(C,[],2);
    Q = reshape(Q,Nt*Ntr,Npt);
    
    % Generate interpolation matrix: pts -> modules
    % W = interpmatrix(C,pts')'; %[Npt,Nm]
    W = interpmatrix(C,pts','extrap','nearest')';
    
    M = zeros([Nm,Nt,Ntr],'single');
    S = zeros([Nm,Nt,Ntr],'single');
    
    % Divide time-steps into chunks, to avoid memory issues
    Nchunks = ceil(Nt*Ntr*Nm*Ns/maxarraysize()*8);
    ChunkSize = ceil(Nt*Ntr/Nchunks);
    for k = 1:Nchunks
        idx = false(Nt*Ntr,1);
        idx((k-1)*ChunkSize + 1 : min(k*ChunkSize,Nt*Ntr)) = true;
        
        % Interpolate at cell centers
        m = reshape(double(Q(idx,:))*W,nnz(idx),Nm,Ns);
        M(:,idx) = mean(m,3)';
        S(:,idx) = std(m,1,3)'; % variation accross cells
    end
    M = shiftdim(M,1);
    S = shiftdim(S,1);
    
%     % Get module centers
%     C(:,1) = arrayfun(@(e) mean(e.border.x),geom.elements);
%     C(:,2) = arrayfun(@(e) mean(e.border.y),geom.elements);
%     Nm = size(C,1);
%     
%     [Nt,Ntr,Npt] = size(Q);
%     assert(Npt == size(pts,2),'Dpoa size does not match node-point coordinates');
%     
%     % Generate interpolation matrix: pts -> modules
%     % W = interpmatrix(C,pts')'; %[Npt,Nm]
%     W = interpmatrix(C,pts','extrap','nearest')';
%     
%     % Interpolate at each module's center
%     Q = reshape(Q,Nt*Ntr,Npt);
%     M = single(reshape(double(Q)*W,Nt,Ntr,Nm));
%     S = 0;
end

function key = checkmodType(modType)
% Parse known aliases into pvl_

    % TODO: Unify material key names (also used in CHECKODM)
    notquite = false;
    switch lower(modType)
    case {'c-si','csi','monosi'}, key = 'monosi';
    case {'p-si','psi','si-poly','polysi'}, key = 'polysi';
    case {'a-si','asi'}, key = 'asi';
    case {'cdte'}, key = 'cdte';
    case {'cigs'}, key = 'cigs';
    case {'cis'}, key = 'asi'; notquite = true;
    case {'hit'}, key = 'monosi'; notquite = true;
    % case {'gaas','gainp','gainas','ge','gain',..}, key = '';
    otherwise, key = '';
    end
    
    if isempty(key) && ~isempty(modType)
        warning('Unknown material: %s',modType);
    end
    if notquite
        warning('Using approx. spectral correction (%s) for %s',key,modType);
    end
end

