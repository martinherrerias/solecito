function varargout = infinite_rows(P,L,tilt,m,varargin)
% [SFA,F_x2s,eqSA,F] = INFINITE_ROWS(P,H,TILT,[M]) - Infinite line (2D) isotropic-diffuse view
%    factors to sky F_x2s, and view-factor matrix F for M sections (rows) of tables with width H,
%    pitch P, and given TILT, following Tschopp et al. [3],[1],[2].
%
%      Dx = F_x2s.*Dh;
%      Ix = f_x.*In;
%      GTI = (Dx + Ix)*inv(I - RHO.*F);
%
%    Use 'nr',A,'ng',B,.. to set the number of ground and module-back reflecting segments. If not
%    specified they will be set to max(6,ceil(M/4)).
%
% [Ge,BshF,DTI,BTI0,GTIm,F_x2s,f_x,Fi,Dh,In,CSn,f_cs] = INFINITE_ROWS(P,H,TILT,M,AZ,DATA,...) 
%   Apply calculated view factors to irradiance and solar position DATA. If differs from [3] 
%   in the treatment of circumsolar fraction (set 'model','hay','CSR',0,.. for the original):
%
%      Dx = F_x2s.*Dh;
%      Ix = f_x.*(In-CSn) + f_cs.*CSn;
%      Fi = inv(I - ALBEDO.*F);
%      Ge = (Dx + Ix)*Fi(1:m,:)';
%
% INPUT:
%   P - pitch (distance among rows, m) - (H + D, in [1]).
%   H - table width (along module surface)
%   TILT - (degrees from horizontal)
%   AZ - (degrees) convention must match DATA.sunaz
%
%   DATA.sunel - [N,1] (degrees) apparent solar elevation
%   DATA.sunaz - [N,1] (degrees) solar azimuth
%   DATA.DHI, BNI, ENI... Diffuse horizontal, Direct-normal, and Extraterrestrial-Normal Irradiance
%
% 'NAME',VALUE options:
%
%   'albedo', 0.2 - albedo (isotropic, lambertian)
%   'fIAM', '' - Cosine response model, parsed through CHECKIAM
%   'CSR', 25 - Circumsolar radius (degrees). Special cases are 0 (dot), and Inf (isotropic).
%   'backrho', 0.9 - back module reflectance
%   'frontrho', 0.1 - front module reflectance
%   'model', 'hay' - 'hay' or 'perez', affects the estimation of circumsolar-fraction
%   'nr', ceil(sqrt(M)) - number of back segments
%   'ng', ceil(sqrt(M)) - number of ground segments
%   'clearance', 0.5 - distance from ground to lower edge of modules
%
%   '-GTI' or 'GTI',true - requires a field DATA.GTI that contains non-shaded tilted irradiance
%       data for the given TILT and AZ. It then applies a linear correction factor to the 
%       calculated Dh (and if necessary In, CSn) so that the transposition matches DATA.GTI.
%
% OUTPUT:
%   SFA - shade-free angle (degrees)
%   F_x2s - [1,Q] sky view factors for front, rear, and ground, F_x2s = [F_cs, F_gs, F_rs]
%           Q = M + Nr + Ng is the total number of reflecting surfaces used in the calculation.
%   eqSA - [Q,2] equivalent shade angles for front, rear, and ground segments. eqSA(:,1)
%       is the angle cast by the front row, and eqSA(:,2) the angle cast by the back (self) row,
%       such that they result point-wise equivalent sky-view factors.
%   F - [Q,Q] full view factor matrix.
%
%   Ge - [N,M] effective tilted irradiance (with linear shading & IAM correction)
%   BshF - [N,M] Beam-shading fraction (0 = not shaded, 1 = fully shaded)
%   DTI - [N,M] Linearly shaded diffuse tilted irradiance (IAM corrected)
%   BTI0 - [N,1] Non-shaded direct tilted irradiance (IAM corrected)
%   GTIm - [N,1] unshaded, uncorrected transposition model output
%   f_x - [N,Q] incidence (and shading) factor for direct irradiance over each surface
%   Fi - [Q,Q] Fi = inv(I - ALBEDO.*F)
%   Dh - [N,1] isotropic diffuse component
%   In - [N,1] direct irradiance component
%   CSn - [N,1] circumsolar normal component
%   f_cs - [N,Q] incidence (and shading) factor for circumsolar irradiance.
%
% [1] Appelbaum, J., 2016. Current mismatch in PV panels resulting from different locations of 
% cells in the panel. Solar Energy 126, 264–275. https://doi.org/10.1016/j.solener.2016.01.013
%
% [2] Appelbaum, J., Massalha, Y., Aronescu, A., 2019. Corrections to anisotropic diffuse 
% radiation model. Solar Energy 193, 523–528. https://doi.org/10.1016/j.solener.2019.09.090
%
% [3] Tschopp, D., Jensen, A.R., Dragsted, J., Ohnewein, P., Furbo, S., 2022. Measurement and
% modeling of diffuse irradiance masking on tilted planes for solar engineering applications. 
% Solar Energy 231, 365–378. https://doi.org/10.1016/j.solener.2021.10.083
%
% [4] Varga, N., Mayer, M.J., 2021. Model-based analysis of shading losses in ground-mounted 
% photovoltaic power plants. Solar Energy 216, 428–438. https://doi.org/10.1016/j.solener.2021.01.047
%
% [5] Han, K., Feng, Y.T., Owen, D.R.J., 1970. An Accurate Algorithm for Evaluating Radiative Heat 
% Transfer in a Randomly Packed Bed. CMES 49, 143–162. https://doi.org/10.3970/cmes.2009.049.143
%
% TODO: account for circumsolar contributions from front & back rows.
% TODO: bifacial case.
%
% See also: MISMATCH, PV_SIMPLE

    narginchk(3,Inf);
    if nargin < 4, m = 1; end
    
    opt.data = [];
    opt.surfaz = [];
    opt.albedo = 0.2;
    opt.fIAM = '';
    opt.CSR = 25;
    opt.backrho = 0.9;
    opt.frontrho = 0.1;
    opt.model = 'hay';
    opt.nr = ceil(sqrt(m));
    opt.ng = ceil(sqrt(m));
    opt.clearance = 0.5;
    opt = parseoptions(varargin,{'-GTI'},opt,'dealrest',2);

    validateattributes(P,'numeric',{'scalar','positive','real'});
    validateattributes(L,'numeric',{'scalar','positive','real'});
    validateattributes(tilt,'numeric',{'scalar','nonnegative','real'});
    validateattributes(m,'numeric',{'scalar','positive','real','integer'});
    opt.model = parselist(opt.model,{'hay','perez'}');
        
    [SFA,F_x2s,eqSA,F] = viewfactors(L,P,tilt,opt.clearance,m,opt.nr,opt.ng);

    if isempty(opt.data)
        varargout = {SFA,F_x2s,eqSA,F};
        return;
    end
       
    validateattributes(opt.surfaz,'numeric',{'scalar','nonnegative','real'});
    validateattributes(opt.data,{'struct','table','timetable','MeteoData'},{'nonempty'});
    
    data = opt.data;
    surfaz = opt.surfaz;
    
    switch class(data)
    case 'table', data = table2struct(data);
    case 'timetable', data = timetable2struct(data);
    case 'MeteoData'
        [data,B] = bestimate(data,'fillgaps',false);
        data = completestruct(B,timetable2struct(data.data));
    end
    parsestruct(data,{'BNI','DHI','ENI','sunel','sunaz'},'real','vector','equalsize');
    
    if opt.GTI
        parsestruct(data,{'GTI'},'finite','real','vector','numel',numel(data.BNI));
    end
    
    [f_x,CIA,BshF] = incidence_factors(data.sunaz,data.sunel,surfaz,L,P,tilt,opt.clearance,m,opt.nr,opt.ng);
    
    switch opt.model
    case 'hay', [Dh,In,Rh] = hay_components(data.DHI,data.BNI,data.ENI,data.sunel,opt.albedo);
    case 'perez', [Dh,In,Rh] = perez_components(data.DHI,data.BNI,data.ENI,data.sunel,opt.albedo,tilt);
    otherwise, error('you should not be here');
    end
    
    BTI0 = data.BNI.*max(0,CIA); % not shaded!
    
    % Treatment of circumsolar
    CSn = In - data.BNI;
    if any(CSn < 0)
       CSn = max(0,CSn);
    end
    
    validateattributes(opt.CSR,'numeric',{'vector','nonnegative','real'});
    assert(isscalar(opt.CSR) || numel(opt.CSR) == numel(BTI0),'Unexpected CSR size');
    if all(opt.CSR == 0)
    % Treat as beam shading
        f_cs = f_x;
    elseif all(opt.CSR == Inf)
    % Mix with isotropic
        f_cs = F_x2s.*max(0,CIA);
    else
    % Apply geometric shading factors
        [f_cs,fg0] = circumsolar_factors(opt.CSR,eqSA,data.sunel,data.sunaz,surfaz,tilt,m,opt.nr,opt.ng);
        
        % Renormalize to ensure that CS component matches model for horizontal non-shaded plane
        switch opt.model
        case 'hay', CSn = CSn.*max(0,sind(data.sunel))./fg0;
        case 'perez', CSn = CSn.*max(0.087,sind(data.sunel))./fg0;
        end
        In = CSn + data.BNI;
    end

    % Reference (non-shaded model)
    cb = cosd(tilt);
    GTIm = Dh.*(1+cb)/2 + Rh.*(1-cb)/2 + In.*max(0,CIA);
    
    if opt.GTI
        % Apply correction to diffuse components to match provided GTI
        [Dh,In] = remove_bias(Dh,In,data.BNI,GTIm,data.GTI);
    end

    rho = [repmat(opt.frontrho,1,m) repmat(opt.albedo,1,opt.ng) repmat(opt.backrho,1,opt.nr)];
    Fi = eye(m+opt.ng+opt.nr) - rho.*F;

    Dx = F_x2s.*Dh;
    Ix = f_x.*(In-CSn) + f_cs.*CSn;
    Ge = (Dx + Ix)/Fi';
    
    Ge = Ge(:,1:m);
    Fi = inv(Fi);

    DTI = Ge - BTI0.*(1-BshF);

    opt.fIAM = checkIAM(opt.fIAM);
    
    if ~all(opt.fIAM(0:89) == 1)
        
        IAM = opt.fIAM(acosd(CIA));
        
        % Equivalent tilt angle that would generate the same view factor
        eq_tilt_iso = acosd(2*F_x2s(1:m) - 1);

        [iso_IAM,gnd_IAM] = diffuseIAM(opt.fIAM,eq_tilt_iso);
        
        if all(opt.CSR == 0)
            cs_IAM = IAM;
        elseif all(opt.CSR == Inf)
            cs_IAM = iso_IAM;
        else
            if isscalar(opt.CSR)
                cs_IAM = circumsolarIAM(opt.fIAM,acosd(CIA),opt.CSR);
            else
                cs_IAM = circumsolarIAM_wrapper(opt.fIAM,acosd(CIA),opt.CSR);
            end
        end
        CS = f_cs(:,1:m).*CSn;
        assert(~any(CS < 0,'all') && ~any(CS > DTI,'all'));

        BTI0 = BTI0.*IAM;
        
        DTI = (Dx(:,1:m) - CS).*iso_IAM + CS.*cs_IAM + (DTI - Dx(:,1:m)).*gnd_IAM;
        Ge = BTI0.*(1-BshF) + DTI;
    end

    varargout = {Ge,BshF,DTI,BTI0,GTIm,F_x2s,f_x,Fi,Dh,In,CSn,f_cs};
end

function [SFA,F_x2s,eqSA,F] = viewfactors(L,P,tilt,h,m,r,g)
% View factors, from Hottel rule, as in Tschopp et al. for m rows (1 rear and ground section!)
% F_i2j = factor from i to j, for c = front, s = sky, g = ground, r = rear
        
    sb = sind(tilt);
    cb = cosd(tilt);
    
    D = P - L*cb; % hor. gap between rows
    assert(D > 0,'Overlapping tables!');

    % pts. along table surface, from bottom edge
    pc = linspace(0,L,m+1)'.*[cb,sb];
    
    % pts. along rear, bottom to top
    pr = linspace(L,0,r+1)'.*[cb,sb] - [P,0];
    
    % pts. along ground
    pg = [linspace(-P,0,g+1)',-h*ones(g+1,1)];
    
    self = [pc(1,:);pc(m+1,:)];
    next = self - [P,0];
    sky = [pc(m+1,:);pr(1,:)];
             
    F_c2s = hottelcrossrule(pc,sky);
    F_r2s = hottelcrossrule(pr,sky);
    
    % Used to dettach obstacles from segments, see note ($) in hottelcrossrule.
    SHIFT = [1e-9,0];

    F_g2s = hottelcrossrule(pg,sky, self+SHIFT) + ...
            hottelcrossrule(pg,sky - [P,0], next+SHIFT) + ...
            hottelcrossrule(pg,sky + [P,0], self-SHIFT);
         
    F_c2g = hottelcrossrule(pc,pg,[0,0; 0,-h]+SHIFT) + ...
            hottelcrossrule(pc,pg - [P,0],next);

    F_x2s = [F_c2s; F_g2s; F_r2s]';

    F_g2r = hottelcrossrule(pg,pr) + ...
            hottelcrossrule(pg + [P,0],pr,self) + ...
            hottelcrossrule(pg - [P,0],pr,[-P,0;-P,L] - SHIFT);
        
    F_c2r = hottelcrossrule(pc,pr);

    F = zeros(m+g+r);
    F(1:m,m+1:m+g) = F_c2g;
    F(1:m,m+g+1:end) = F_c2r;
    F(m+1:m+g,1:m) = F_c2g'*(L/m)./(P/g);
    F(m+1:m+g,m+g+1:end) = F_g2r;
    F(m+g+1:end,1:m) = F_c2r'*r/m;
    F(m+g+1:end,m+1:m+g) = F_g2r'*(P/g)/(L/r);

    SFA = atan2d(L*sb - pc(:,2)', pc(:,1)' + D); % shade-free angle / "obscuring angle" [1]
    SFA(end) = [];

    eqSA = equivalentshadeangles(L,P,tilt,F_x2s,h,m,r,g);
    
    function F = hottelcrossrule(XY1,XY2,obs)
    % Hottel cross-string view factor, according to Han et al. 1970. From points XY1 ti XY2,
    % potentially obstructed by obstacle OBS.

        daB = @(a,B) hypot(a(1)-B(:,1),a(2)-B(:,2));
        
        if nargin < 3 || isempty(obs)
            d = pdist2(XY1,XY2,daB);
        else
            d = pdist2(XY1,XY2,@distfun);
        end
        s2 = size(d,2)-1;
        s1 = size(d,1)-1;
        
        dL = vecnorm(XY1(1:s1,:) - XY1(2:s1+1,:),2,2);
        
        F = zeros(s1,s2);
        for j = 1:s1
            F(j,:) = 0.5*(d(j,1:s2) + d(j+1,2:s2+1) - d(j+1,1:s2) - d(j,2:s2+1))'/dL(j);
        end
        assert(all(F > -1e-6,'all'));
        F = max(0,F);
        
        function d = distfun(a,B)
        % Euclidean distance from A to B, unless the ray A-B intersects segment OBS, then return
        % the distance around OBS(1,:), i.e. A-OBS(1,:)-B
        
            TOL = 1e-12;
            d = zeros(size(B,1),1);
        
            [ix,~,t] = segment_intersection(B,a,obs(1,:),obs(2,:),TOL);
            
            % ($) NOTE: Rays A-B can't start/end at obstacle, otherwise it's not clear whether 
            %     the distance should be measured directly or around OBS.
            assert(~any(abs(t(ix)) < TOL | abs(1-t(ix)) < TOL),'Intersection too close to ray end')
            
            d(~ix) = daB(a,B(~ix,:));
            d(ix) = daB(obs(1,:),B(ix,:)) + daB(a,obs(1,:));
        end
    end
end

function [f_x,CIA,BshF,prjel] = incidence_factors(sunaz,sunel,surfaz,H,P,tilt,h,m,r,g)
        
    if nargin == 0, test(); return; end

    sunvec = sph2cartV(90-sunaz,sunel);
    
    % solar elevation projected to the plane perpendicular to the modules & the ground
    prjel = atan2d(tand(sunel),cosd(sunaz-surfaz)); 
    
    % cosine of incidence angle
    CIA = sunvec*sph2cartV(90-surfaz,90-tilt)';

    sb = sind(tilt);
    cb = cosd(tilt);

    % beam-shading fraction, per table and per row
    shf = P/H*(sb./tand(tilt + prjel) - cb).*(CIA > 0) + 1;
    BshF = rowshadingfactors(shf,m);

    fc = max(0,CIA).*(1-BshF);

    fr = P/H*(cb - sb./tand(tilt + prjel)).*(CIA < 0) + 1;
    fr = (1 - rowshadingfactors(fr,r)).*max(0,-CIA);

    shade_start = mod(1+h./tand(prjel)/P,1);
    shade_stops = mod((H*cb + (H*sb+h)./tand(prjel))/P,1);
    flipped = shade_start > shade_stops;
                  
    shg = (1 - rowshadingfactors(shade_start,g)) + ...
          rowshadingfactors(shade_stops,g) - 1*(~flipped);
                  
    allshaded = prjel < atan2d(H*sb,P-H*cb) | prjel > atan2d(H*sb,-P-H*cb);
    shg(allshaded,:) = 1;
    
    fg = (1-shg).*(max(0,sunvec(:,3)));

    f_x = [fc,fg,fr];
end

function BshF = rowshadingfactors(shf,m)
    
    if isempty(shf), BshF = zeros(0,m); return; end
    
    shf = double(min(1,max(0,shf)));
    allshaded = (shf*m >= (1:m));
    partshaded = ~allshaded & (shf*m + 1 > (1:m));
    BshF = allshaded + partshaded.*rem(shf*m,1);
    
    assert(all(abs(sum(BshF,2) - shf*m) < eps(m)));
end

function [Dh,In,Rh] = hay_components(DHI,BNI,ENI,sunel,albedo)
        
    cz = sind(sunel);

    Ai = BNI./ENI;
    Dh = DHI.*(1-Ai);
    In = BNI + Ai.*DHI./cz;

    GHI = DHI + max(0,BNI.*cz);
    Rh = GHI.*albedo;
end

function [Dh,In,Rh] = perez_components(DHI,BNI,ENI,sunel,albedo,tilt)
        
    cb = cosd(tilt);

    [F1,F2] = pvlmod_perezcoeffs(BNI,DHI,ENI,sunel);
        
    cz = sind(sunel);
    b = max(0.087,cz); 
    
    Fiso = (1 + cb)/2; 
    Fhb = sind(tilt);

    % Merge isotropic with HBB/Fiso, so that D = Dh Fiso = DHI(1 - F1) Fiso + F2 Fhb 
    Dh = DHI.*(1 - F1 + F2.*Fhb./Fiso);
    In = BNI + F1.*DHI./b;

    GHI = max(0,DHI + BNI.*cz);
    Rh = GHI.*albedo;
end

function [f_cs,fg0] = circumsolar_factors(CSR,eqSA,sunel,sunaz,surfaz,tilt,m,r,g)
% Apply geometric circumsolar shading factors, similar to Varga & Mayer, but accounting for the
% compound effect of several shading planes (front row, back row, horizon), to avoid double
% penalization when the same part of the circumsolar disk is shaded.
%
% See also: INFINITE_ROWS
    
    validateattributes(CSR,'numeric',{'positive','<',90});
    assert(size(eqSA,1) == m+r+g);

    if isscalar(CSR)
        [~,cs_vf] = circumsolarIAM(1,[],CSR);
    else
        [~,cs_vf] = circumsolarIAM_wrapper(1,[],CSR);
    end
    
    sunvec = sph2cartV(90-sunaz,sunel);
    vanishpt = sph2cartV(-surfaz,0)';
    
    surfnormal = sph2cartV(90-surfaz,90-tilt)';

    % Sine of a pseudo-hour-angle (90° if the sun sets at the vanishing point)
    % modulates the approx. angle between shadows cast by two planes that intersect at vanishpt
    sinw = sunvec*vanishpt;
        
    % solar elevation projected to the plane perpendicular to the modules & the ground
    prjel = atan2d(tand(sunel),cosd(sunaz-surfaz)); 
    
    % cosine of incidence angle
    CIA = sunvec*surfnormal;
        
    % Approx angle between shading semiplanes (front row @ eqSA, and ground) 
    phi = eqSA(1:m,1)'.*sinw;
    
    % Shading due to front row and ground
    fc = circumsolarshading(prjel-eqSA(1:m,1)',sunel,phi,CSR,'-deg');

    % Incidence + self-shading
    fc = cs_vf(acosd(CIA)).*fc;

    % Shading of back surface
    phi = eqSA(m+g+1:end,2)'.*sinw;
    fr = circumsolarshading(eqSA(m+g+1:end,2)'-prjel,sunel,phi,CSR,'-deg');
    fr = cs_vf(acosd(-CIA)).*fr;

    % Shading of ground
    phi = permute(eqSA(m+1:m+g,:),[3 1 2]).*sinw;
    
    fg0 = cs_vf(90-sunel); % incidence + self
    
    s0 = circumsolarshading(sunel,CSR,'-deg');
    fg = circumsolarshading(prjel - eqSA(m+1:m+g,1)',sunel,phi(:,:,1),CSR,'-deg')./s0.* ... % front row
         circumsolarshading(eqSA(m+1:m+g,2)'-prjel,sunel,phi(:,:,2),CSR,'-deg')./s0;        % back row
    fg = fg0.*fg;                                              

    f_cs = [fc,fg,fr];
end

function [iam,cs_vf] = circumsolarIAM_wrapper(fIAM,IA,CSR)
    CSR = round(CSR);
    [r,~,ia] = unique(CSR);
    inj = arrayfun(@(j) (ia == j),1:numel(r),'unif',0);
    if isempty(IA)
        for j = numel(r):-1:1
            [iaf{j},csf{j}] = circumsolarIAM(fIAM,[],r(j));
        end
        cs_vf = @(x) vectorfun(csf,x);
        iam = @(x) vectorfun(iaf,x);
    else
        cs_vf = zeros(size(IA));
        iam = zeros(size(IA));
        for j = numel(r):-1:1
            [iam(inj{j}),cs_vf(inj{j})] = circumsolarIAM(fIAM,IA(inj{j}),r(j));
        end
    end

    function vf = vectorfun(f,x)
        vf = zeros(size(x));
        for k = 1:numel(f)
            vf(inj{k}) = f{k}(x(inj{k}));
        end
    end
end
    
function eqSA = equivalentshadeangles(H,P,tilt,F_x2s,h,m,r,g)
% Equivalent shade angles for each surface, corresponding to points that have the same view-factor
% as the integrated F_x2s, calculated for the whole surface (by VIEWFACTORS, above).
%
% TODO: equivalent shade angles for ground segments include only the center 'gap' of sky on top of
% the ground. For segments close & below the modules, contribution from neighboring sky gaps can
% be more important (See T20221020_ground_view_factor_infinite_line.m)
% eqSA should contain at least 6 equivalent angles for these cases.

    assert(m + r + g == numel(F_x2s));

    eqSA = zeros(m+r+g,2);
    eqSA(1:m,1) = acosd(2*F_x2s(1:m) - 1) - tilt;
    eqSA(1:m,2) = 180 - tilt;

    hy = H*sind(tilt)+h;
    hx = H*cosd(tilt);
    
    % Ground: find a point within each interval that matches Fg2s, see note above
    VF = @(x) 0.5*(min((P+hx-x)./hypot(P+hx-x,hy),(P-x)./hypot(P-x,h)) - (hx-x)./hypot(hx-x,hy)) + ...
              0.5*max(0,(2*P+hx-x)./hypot(2*P+hx-x,hy) - max((P+hx-x)./hypot(P+hx-x,hy),(P-x)./hypot(P-x,h))) + ...
              0.5*max(0,(x+P-hx)./hypot(x+P-hx,hy) - x./hypot(x,h)); 

    gg = linspace(0,P,g+1);
    [x,bounded] = bisection(@(x) VF(x) - F_x2s(m+1:m+g),gg(1:g),gg(2:g+1));
    x(~bounded) = 0.5*(gg([~bounded,false]) + gg([false,~bounded]));
    
    eqSA(m+1:m+g,1) = atan2d(hy,x-hx);
    eqSA(m+1:m+g,2) = 180 - atan2d(hy,P+hx-x);
    
    % Rear
    eqSA(m+g+1:end,1) = 180 - tilt; 
    eqSA(m+g+1:end,2) = 360 - acosd(2*F_x2s(m+g+1:end) - 1) - tilt;
end

function [Dh,In] = remove_bias(Dh,In,GTI_model,GTI_ref)
% Apply correction to diffuse components to match provided GTI

    weird = 1 - In./GTI_ref < 0.05;
    k = (GTI_ref - In)./(GTI_model - In);
    k(weird) = GTI_ref(weird)./GTI_model(weird);
    Dh = k.*Dh;
    In(weird) = k(weird).*In(weird);
end