function varargout = system_designer(MD,SP,Loc,ModIV,celltemp,Inverter,GCR,mismatch)
% [TILT,AZ,MPS,NS] = SYSTEM_DESIGNER(MD,SP,Loc,ModIV,celltemp,Inverter,GCR)
%   Estimate optimum array TILT and azimuth (AZ), modules per string (MPS) and strings per inverter
%   (NS) for an area ocupation rate GCR (default 0.4), meteo-data MD, solar-position SP, location
%   LOC, and module, cell-temperature, and Inverter models MODIV, CELLTEMP, INVERTER.
%
%   Simplifying assumptions are infinite-row shading (diffuse isotropic), and linear shading effect

    if nargin < 7
        warning('Assuming Ground-Cover-Ratio of 0.4');
        GCR = 0.4;
    end
    if nargin < 8, mismatch = 'h'; end
    
    eta_STC =ModIV.getMPP(1000,25)./(1000*ModIV.area); % (Wp/m²)/(1000 W/m²) = Wp/kW
    fprintf('Required area: %0.1f m²/kWp\n',1/(eta_STC*GCR));

    Trck.type = '0a';
    % Trck.azimuth = 0;
    %Trck.tilt = EL;
    Trck.bankangle = 0;
    Trck.tracklimits = [-45,45];
    Trck.backtracking = true;
    Trck.groundcoverratio = GCR;

    % Perez et al. anisotropy indices for Circumsolar and Horizon-Brightening
    % [F1,F2] = pvlmod_perezcoeffs(MD.BNI,MD.DHI,MD.ENI,SP.El);
    % s = sph2cartV(SP.Az,SP.El);
    % b = max(0.087,s(:,3));  % b = max(0.087,k's)
    
    % [t,a] = meshgrid(0:5:90,-45:5:45);
    % GUIfigure('system_designer_tilt','Optimum Orientation');
    % clf(); hold on; grid on;
    % g = arrayfun(@(t,a) mean(tilted(a,t),'omitnan'),t,a);
    % contourf(t,a,g); colorbar;
    % xlabel('Tilt');
    % ylabel('Azimuth');
    
    x0 = [0,min(20,Loc.latitude)];
    f = @(x) -mean(tilted(x(1),x(2)),'omitnan');
    opt = optimset();
    opt.MaxFunEvals = 2000;
    x0 = fminsearch(f,x0,opt);
    fprintf('Optimum orientation: %0.1f° tilt, %0.1f° azimuth\n',x0(2),x0(1));
    
    GUIfigure('system_designer_tilt','Optimum Orientation');
    clf(); hold on; grid on;
    t = unique(max(0,min(90,(-10:2:10) + round(x0(2)))));
    a = (-10:2:10) + round(x0(1));
    [t,a] = meshgrid(t,a);
    g = arrayfun(@(t,a) mean(tilted(a,t),'omitnan'),t,a);
    contourf(t,a,g); colorbar;
    xlabel('Tilt');
    ylabel('Azimuth');
    
    [~,GTI] = tilted(round(x0(1)),round(x0(2))); % assigns tilt = x0(2), az = x0(1)
    plot(Trck.tilt,Trck.azimuth,'r+');

    if isempty(Inverter), varargout = {Trck.tilt,Trck.azimuth}; return; end
        
    % Optimize modules-per-string (mps) and strings per inverter (ns)
    Tc = celltemp(GTI,MD.Ta,MD.vw);
    [P0,V0] = ModIV.getMPP(GTI,Tc);
    mps = [floor(Inverter.MPPTLow/prctile(V0,25)),ceil(Inverter.MPPTHi/prctile(V0,90))];
    ns = Inverter.Pdc0*[0.5;1.2]./(mps*prctile(P0,95));
    ns = [floor(min(ns(:))),ceil(max(ns(:)))];
    ns = max(1,ns);
    
    [mps,ns] = meshgrid(mps(1):mps(end),ns(1):ns(end));
    Pac = arrayfun(@(m,s) mean(ONDeval(Inverter,P0*m*s,V0*m,MD.Ta),'omitnan'),mps,ns);
    eff = Pac./(mean(P0,'omitnan').*mps.*ns);
    [~,idx] = max(eff(:));
    
    GUIfigure('system_designer_elec','String configuration');
    clf(); hold on; grid on;
    lvl = sort(1-logspace(log10(1-max(eff(:))),log10(1-median(eff(:))),10));
    contourf(mps,ns,eff*100,lvl*100);
    xlabel('Modules per string');
    ylabel('Strings per MPPT');
    plot(mps(:),ns(:),'ko');
    plot(mps(idx),ns(idx),'ro','markerfacecolor','r');
    colorbar('ticks',lvl*100);
    title('Inverter efficiency');
    
    varargout = {Trck.tilt,Trck.azimuth,mps(idx),ns(idx)};
    
    % Mod = importpolygonsfile('*.mpoly');
    % [w,h] = rectangleproperties(Mod.geom.border);
    % samplesystem('0a',[w,h],[n, m, N, M],[sn,sm,sN,sM],'az',Trck.azimuth,'tilt',Trck.tilt)
    
    function [P0,GTI] = tilted(a,t)
        
        Trck.azimuth = a;
        Trck.tilt = t;

        u = mountrotations(Trck,SP.Az,SP.El);   % [3,3,Nm,Nt]
        St = squeeze(atan2d(hypot(u(1,3,:,:),u(2,3,:,:)),u(3,3,:,:)));
        Sa = squeeze(atan2d(u(2,3,:,:),u(1,3,:,:)));
        
        % Beam Shading Factor, infinite rows
        F_b = (1/GCR)./(cosd(St) + sind(St).*tand(90-SP.El).*cosd(Sa - SP.Az));
        F_b = max(0,min(1,F_b));
        F_b(SP.El <= 0) = 0;
        
        switch lower(mismatch)
            case {'h','linear'}
            case {'v','worst'}, F_b(F_b < 71/72) = 0;
            otherwise
                error('Unknown mismatch condition');
        end
        
        % Isotropic shading, infinite rows
        shfa = atand(GCR.*sind(St)./(1-GCR.*cosd(St))); % shade-free angle
        F_d = 0.5*(1+cosd(St + shfa));
        
        [~,ISO,CS,HB,ALB,BTI] = pvlmod_perez(St,Sa,MD.DHI,MD.BNI,MD.ENI,SP.El,SP.Az,MD.albedo);
        GTI = (CS + BTI).*F_b + (HB + ALB + ISO).*F_d;
        
        Tc = celltemp(GTI,MD.Ta,MD.vw);
        P0 = ModIV.getMPP(GTI,Tc,'interp');
    end
end