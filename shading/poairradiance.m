function [POA,ShR] = poairradiance(MD,SP,Trck,horizon,ShR,varargin)
% POA = POAIRRADIANCE(MD,SP,TRCK,HOR,[SHR],..) - Calculate Plane-Of-Array irradiances
%   for the Nt time-steps of of Meteo-Data MD and Solar-Position SP; and for the Nu analysed
%   mounts in TRCK, considering far-horizon HOR, and (optionally) shading-results SHR.
%
%   POAIRRADIANCE combines the view-factors included in SHR with a radiance-distribution estimated
%   by DIFFUSECOMPONENTS
%
% INPUT:
%   MD - structure with Nt-vector fields BNI (Beam Normal-), DHI (Diffuse Horizontal-), GHI 
%		(Global	Horizontal-), ENI (Extraterrestrial Irradiance), albedo, and optionally AMr 
%       (relative air mass), and custom Perez coefficients F1, F2. All Nt-vectors.
%   SP - structure with fields Az (azimuth, Equator = 90° convention) and El (elevation), both 
%		Nt-vectors in degrees.
%   TRCK - mounts structure, the fields required depend on TRCK.type (see MOUNTROTATIONS)
%       The optional field TRCK.analysedtrackers (default 1:end) can reduce/resort the mounts
%       under analysis.
%   HOR - (spherical) POLYGON3D object, far-horizon profile (common to sensors and mounts).
%   SHR - (optional) SHADINGRESULTS object (step-by-step or for interpolation)
%
%  ..,'model',MODEL - for MODEL in {'perez','haydavies','igawa'}, determines 
%       whether to use the Perez et al. 1990, Hay & Davies 1980, or Igawa et al. 2004 models for 
%       irradiance transposition. The default is SimOptions.diffusemodel.
%
%       To specify a non-standard set of coefficients for the perez model, use the syntax 
%       'perez - SETNAME' where SETNAME can be any of the sub-models in PVLMOD_PEREZ.
%       Alternatively, if fields F1 and F2 are passed along with MD, they will be used as custom
%       anisotropy indices for circumsolar-fraction and horizon-brightening.
%
%   ..,'geom',RG - Provide a SHADINGREGIONS object specifying the geometrical framework for
%       irradiance decomposition. Currently, when using Perez or Hay-Davies models, RG must reflect
%       the corresponding model, i.e. using SHADINGREGIONS('perez'). 
%       Igawa takes any geometry. The default is SHR.worldgeom, or SimOptions.skyregions.
%
%   ..,'fiam',IAM - Override IAM function. Default is SHR.info.options.IAM, or just checkIAM()
%
% OUTPUT: in the following, Nu,Nt,Nm,Np have the meaning used in SHADINGRESULTS
%
%   POA.Bpoa0: Nt·Nu* array of tilted beam-irradiance values, considering far-horizon only, and
%       perfect absortivity (i.e. IAM = 1)
%   POA.Dpoa0: Nt·Nu* array of tilted diffuse-irradiance values (id.)
%
%   POA.Bpoa_IAM: Nt·Nu* array of tilted beam-irradiance values, considering far-horizon and IAM.
%   POA.Dpoa_IAM: Nt·Nu* array of tilted diffuse-irradiance values (id.)
%   POA.CS: Nt·Nu*·Np* array of tilted (shaded, IAM-adjusted) circumsolar-irradiance values. 
%       If SHR is not available, or it doesn't contain circumsolar regions, then POA.CS will be 
%       an Nt·Nu* estimate (linear shading, Perez model).
%
%   POA.Dpoa: Nt·Nu·Np array of tilted (shaded) diffuse-irradiance values, using SHR
%   POA.ShF: Nt·Nu·Nm array, expanded SHR.BshF with fractions of shade-coverage for each module 
%       (1..Nm) in every mount (1..Nu), at every timestep. 1 = full shade, 0 = no shade.
%   POA.Nsb: Nt·Nu·Nm array, expanded SHR.Nsb with number of shaded cell-blocks on every module.
%
%   POA.bias: Nt vector of correction factors applied to sky- and circumsolar-components to account
%       for finite sky-region sizes, i.e. to ensure that the product sum of radiance components and
%       view-factors is exact on the horizontal.
%
% See also: PVLMOD_PEREZ, PVLMOD_HAYDAVIES, RADIANCE_IGAWA, SHADINGREGIONS, DIFFUSECOMPONENTS,
%           EFFECTIVEIRRADIANCE, SHADINGRESULTS

    narginchk(3,11);
    if nargin < 4, horizon = []; end
    if nargin < 5, ShR = []; end
    
    opt.model = getSimOption('diffusemodel');
    
    if ~isempty(ShR)
        opt.geom = ShR.worldgeom;
        opt.fiam = ShR.info.options.IAM;
    else
        opt.fiam = [];
        opt.geom = ShadingRegions(getSimOption('skyregions'));
    end
    opt = getpairedoptions(varargin,opt,'restchk');
    opt.fiam = checkIAM(opt.fiam);

    assert(isa(opt.geom,'ShadingRegions'),'Expecting SHADINGREGIONS object');
    assert(ischar(opt.model),'Expecting character-array MODEL');    
    parsestruct(MD,{'BNI','DHI','GHI','ENI','albedo'},'-n','-r','-v','-e');
    parsestruct(SP,{'El','Az'},'-n','-f','-v','size',[numel(MD.GHI),1]);
    parsestruct(Trck,{'centers'},'-n','-r','-m','size',[3,NaN]);

    Nu = size(Trck.centers,2);  
    if isfield(Trck,'analysedtrackers') && ~isequal(Trck.analysedtrackers(:)',1:Nu)
        Trck = filterstructure(Trck,Trck.analysedtrackers,Nu);
        Trck = filterstructure(Trck,Trck.analysedtrackers,Nu,'dim',2);
        Nu = size(Trck.centers,2);
    end
    
    mdl = lower(opt.model);
    if contains(mdl,'iso')
        if contains(mdl,'hard')
            warning('Using isotropic shading on Perez components');
            [POA,ShR] = isotropicPOA(Trck,MD,SP,ShR,'hard');
        else
            warning('Using isotropic+circumsolar shading on Perez components');
            [POA,ShR] = isotropicPOA(Trck,MD,SP,ShR,'soft');
        end
        return;
    end
    
    % Unit Vectors pointing towards the sun (project coords. X = 0, Y = 90)
    s = sph2cartV(90 - SP.Az - Trck.rotation,SP.El)'; %[3,Nt]
    
    if ~isempty(ShR)
    % Parse, [interpolate] and unpack shading results
        assert(ShR.Nu == Nu && isequal(Trck.geom.dims,ShR.mountgeom.dims),...
            'Shading Results do not match Trackers structure');
        
        if ~isequal(SP.Az,ShR.az) || ~isequal(SP.El,ShR.el)
            ShR = interpolate(ShR,SP.Az,SP.El);
        end
        
        % Unpack (and keep only) DwF0, DwF, DshF, BshF, and Nsb
        S = struct(ShR,'-unpack','only',{'DwF0','DwF','DshF','BshF','Nsb'});
    end

    % Resolve the problematic cases MD.BNI > 0 & belowhorizon
    if isempty(horizon)
        belowhorizon = SP.El < 0; % technically -32', but cos(z) will make it 0 anyways
    else
        AzPrj = polyprojector('azim','clip',false);
        AzPrj.angtol = min(AzPrj.angtol,32/60); % resolution < sun-diameter
        prjsp = AzPrj.prj(s);
        prjhor = AzPrj.project(horizon,[0;0],true);
        belowhorizon = ~inpolygon(prjsp(1,:),prjsp(2,:),prjhor.x,prjhor.y)';
    end
    trouble = (MD.BNI > 0) & belowhorizon;
    if any(trouble) 
        MD.BNI(belowhorizon) = 0; 
        MD.DHI(belowhorizon) = MD.GHI(belowhorizon);
        warning(['%d time-steps with sun below horizon have non-zero BNI ',...
                 '(values up to %0.1f W). They will be set to DHI = GHI, BNI = 0.'],...
                 nnz(trouble),max(MD.BNI(trouble)));
    end
    
    % Rotation Matrices for trackers, at all time-steps, according to type: [3,3,Nu*,Nt*] array
    n = mountrotations(Trck,SP.Az,SP.El);
    % keep only z' unit vectors on a [3,Nt,Nu] array
    n = permute(n(:,3,:,:),[1,4,3,2]);
    
    % Cosine of incidence angle cos(IA) = max(s·n,0)
    CIA = shiftdim(max(0,sum(n.*s,1)),1); % [Nt,Nu] 
    
    % Calculate radiance distribution
    [D.sky,D.albedo,D.solar] = diffusecomponents(opt.geom,MD,SP,opt.model); % [Nt,Nc]
    
    prj = polyprojector('ortho','normalize',false); % IAM = 1
    A0 = opt.geom.solidangles;

    % Reference (no shading, IAM = 1)
    [F0.sky,F0.albedo,F0.solar] = viewfactors(opt.geom,prj,[0;0;1],s,horizon); % [Nt,Nc]
    [DHIm,~,GND] = applyviewfactors(F0,D,A0); % modeled DHI
    
    % Bias correction
    POA.bias = (MD.DHI-GND)./(DHIm-GND);     % measured / model (sky- and circumsolar-regions)
    POA.bias(~isfinite(POA.bias) | POA.bias < 0) = 1;
    D.sky = D.sky.*POA.bias;
    D.solar = D.solar.*POA.bias;
    
    % For Debugging
    % [tilt,az] = mountangles(Trck,SP.Az,SP.El);
    % [GTI,ISO,CS,HB,ALB,BTI] = pvlmod_perez(tilt,az,MD.DHI,MD.BNI,MD.ENI,SP.El,SP.Az,MD.albedo);

    % No shading, IAM = 1
    [F.sky,F.albedo,F.solar] = viewfactors(opt.geom,prj,n,s,horizon); % [Nt,Nc,Nu]
    POA.Dpoa0 = applyviewfactors(F,D,A0);
    POA.Bpoa0 = single(MD.BNI.*CIA);
    
    % Estimate circumsolar component according to Perez
    b = max(0.087,s(3,:)');
    F1 = pvlmod_perezcoeffs(MD.BNI,MD.DHI,MD.ENI,SP.El);
    POA.CS = F1.*MD.DHI.*CIA./b; % F1·DHI·a/b

    IAM = opt.fiam(acosd(CIA));
    if isempty(ShR) && all(IAM(CIA > 0) == 1)
       POA.Bpoa_IAM = POA.Bpoa0;
       POA.Dpoa_IAM = POA.Dpoa0;
       POA.CS = single(POA.CS);
       return;
    end
    POA.Bpoa_IAM = single(MD.BNI.*CIA.*IAM);
    POA.CS = single(MD.BNI.*CIA.*IAM);
    
    if isempty(ShR)
    % No shading, IAM projection
        prj = polyprojector('IAM','fiam',opt.fiam,'normalize',false);
        [F.sky,F.albedo,F.solar] = viewfactors(opt.geom,prj,n,s,horizon); % [Nt,Nc,Nu]
        POA.Dpoa_IAM = applyviewfactors(F,D,A0);
    else
    % Shaded & not-shaded, IAM projection
    
        if isfield(S.DwF,'gndbeam') && ~isempty(S.DwF(1).gndbeam)
            rd = min(1,max(0.1,(1-F1).*MD.DHI./MD.GHI));
            D.gndbeam = D.albedo.*(1-rd)./b;
            D.albedo = D.albedo.*rd;
            A0.gndbeam = A0.albedo;
        else
            S.DwF = rmfield(S.DwF,'gndbeam');
            S.DshF = rmfield(S.DshF,'gndbeam');
        end
    
        POA.Dpoa_IAM = applyviewfactors(S.DwF,D,A0);
        
        [POA.Dpoa,~,~,CS] = applyviewfactors(S.DshF,D,A0);
        
        if any(CS > 0,'all')
            POA.CS = CS;
        else
            POA.CS = (1-S.BshF).*POA.CS; 
        end
        
        POA.ShF = S.BshF;
        POA.Nsb = S.Nsb;
    end
end

function [POA,ShR] = isotropicPOA(Trck,MD,SP,ShR,type)
% POA = ISOTROPICPOA(TRCK,MD,SP,SHR,TYPE) - Applies shading factors from SHR onto fully isotropic
%   (TYPE = 'hard') or circumsolar+isotropic (TYPE = 'soft') components calculated directly from
%   the Perez model. Input/Output is meant to be equivalent to POAIRRADIANCE.

    A0 = ShR.worldgeom.solidangles;
    ShR = interpolate(ShR,SP.Az,SP.El);

    % Unpack (and keep only) DwF0, DwF, DshF, BshF, and Nsb
    S = struct(ShR,'-unpack','only',{'DwF0','DwF','DshF','BshF','Nsb'});
    
    % Estimate isotropic static factors as the ratio of a weighted unit-radiance,
    % before and after shading
    dummy = struct('sky',1,'albedo',1,'solar',1);
    [~,Asky0,Agnd0,Acs0] = applyviewfactors(S.DwF,dummy,A0);
    [~,Asky,Agnd,Acs] = applyviewfactors(S.DshF,dummy,A0);
    
    F_iso = Asky./Asky0; F_iso = max(0,F_iso); % remove NaNs
    F_alb = Agnd./Agnd0; F_alb = max(0,F_alb);
    if any(Acs,'all')
        F_cs = Acs./Acs0; F_cs = max(0,F_cs);
    else
        F_cs = 1;
    end

    % Get unshaded Perez components
    [t,a] = mountangles(Trck,SP.Az,SP.El);
    [~,ISO,CS,HB,ALB,BTI] = pvlmod_perez(t,a,MD.DHI,MD.BNI,MD.ENI,SP.El,SP.Az,MD.albedo);

    fIAM = checkIAM(ShR.info.options.IAM);
    [IAM_iso,IAM_alb,IAM_hb] = diffuseIAM(fIAM);
    
    CIA = BTI./MD.BNI;
    CIA(MD.BNI <= 0) = 0;
    IAM_beam = fIAM(acosd(CIA));
    
    POA.Bpoa0 = single(BTI);
    POA.Bpoa_IAM = single(BTI.*IAM_beam);
    POA.Dpoa0 = single(ISO + CS + HB + ALB);
    
    POA.ShF = S.BshF;
    POA.Nsb = S.Nsb;
    
    switch lower(type)
    case 'hard'
        POA.Dpoa_IAM = POA.Dpoa0.*IAM_iso(0);
        POA.Dpoa = single(((ISO + HB + CS).*F_iso + ALB.*F_alb)*IAM_iso(0));
    case 'soft'  
        ISO = ISO.*IAM_iso(t);
        ALB = ALB.*IAM_alb(t);
        HB = HB.*IAM_hb(t);
        CS = CS.*IAM_beam;
        POA.Dpoa_IAM = single(ISO + HB + ALB + CS);
        POA.Dpoa = single((ISO + HB).*F_iso + ALB.*F_alb + CS.*F_cs);
    end
    
    POA.bias = NaN(size(POA.Bpoa0,1),1);
end

function varargout = applyviewfactors(F,D,A)
% [DPOA,C1,C2,..] = APPLYVIEWFACTORS(F,D,A) - Apply viewfactors F (multiplied by solid-angles A)  
%   to diffuse components D. F,D,A are all structures with the same fields f = {C1,C2,..}.
%   e.g. {sky,albedo,solar}
% 
%   D.(f) contains [Nt,Nc(f)] matrices of diffuse-radiances [W/m²sr], for Nt timesteps and Nc(f)
%       components/regions of each 'type' (i.e. field) f.
%   A.(f) contains [1,Nc(f)] vectors of solid-angles (stereoradians) for each component/region
%   F can be a [p,q,..] array of structures with equal-sized [Nt,Nc] fields (e.g. from
%     ShadingResults.DshF), or a scalar structure with [Nt,Nc,p,q,..] fields. Typically [p,q,..]
%     is just [Nu,Np] for Np shading-analysis points on Nu mounts.
%
%   APPLYVIEWFACTORS will sum accross 'components' or 'regions' of the same type (i.e. within the
%   same field to return [Nt,p,q,..] arrays C1,C2, of integrated plane-of-array irradiance 
%   values of each type, plus the global integrated DPOA = C1 + C2 + ..
%
% See also: POAIRRADIANCE, ISOTROPICPOA

    narginchk(3,3);

    FLD = fieldnames(F);
    n = numel(FLD);

    assert(isstruct(F) && isstruct(F) && isstruct(A),'Expecting structures F, D, A');
    parsestruct(F,FLD,'-n','-r','-p','-f');
    parsestruct(A,FLD,'-n','-r','-p','-f');
    parsestruct(D,FLD,'-n','-r','-p');

    if ~isscalar(F)
    % Turn F(j,..).x into F.x(:,:,j,..)
        stack = @(varargin) reshape(cat(ndims(varargin{1})+1,varargin{:}),[size(varargin{1}),size(F)]);
        F = struct2cell(F); % [n,p,q] cell of [Nt,Nc] matrices
        F = arrayfun(@(j) stack(F{j,:}),1:n,'unif',0);
        F = cell2struct(F(:),FLD);
    end

    varargout = cell(n+1,1);
    varargout{1} = 0;
    for j = 1:n
       compatiblesize(F.(FLD{j}),D.(FLD{j}),A.(FLD{j}));
       
       if isempty(F.(FLD{j}).*A.(FLD{j}))
        % Don't expand everything just because F.(f) is [Nt,0,..]
            F.(FLD{j})(2:end,:) = [];
            A.(FLD{j})(2:end,:) = [];
       end
        
       varargout{j+1} = single(sum(D.(FLD{j}).*F.(FLD{j}).*A.(FLD{j}),2));      % [Nt,1,p,q,..]
       varargout{j+1} = permute(varargout{j+1},[1,3:ndims(varargout{j+1}),2]);  % [Nt,p,q,..]
       varargout{1} = varargout{1} + varargout{j+1};
    end
end