function [Dsky,Dgnd,Dcs] = diffusecomponents(SR,MD,SP,model)
% [DSKY,DGND,DCS] = DIFFUSECOMPONENTS(SR,MD,SP,MODEL) - Calculate average diffuse radiance for
%   each region of an arbitrary unit-sphere polygonal tesellation SR (SHADINGREGIONS object),
%   given meteo-data MD solar-positions SP, and using an irradiance-distribution/transposition  
%   MODEL {'haydavies','perez[:dataset]','igawa'}.
%
%   DIFFUSECOMPONENTS currently uses GHI,DHI, and BNI (along with AMr and CS_GHI for Igawa model);
%       considering a flat horizon, and making no use of GTI information. 
%
%   FUTURE: Apply inverse-transposition model using several GTI channels and sensor-shading-data.
%
% INPUT:
%   SR - SHADINGREGIONS object specifying the geometrical framework for shading. Currently, when
%       using Perez or Hay-Davies models, SR must reflect the corresponding model, e.g. using
%       SHADINGREGIONS('perez'). Igawa takes any geometry, the default is SimOptions.skyregions.
%
%   MD - structure with Nt-vector fields BNI (Beam Normal-), DHI (Diffuse Horizontal-), GHI 
%		(Global	Horizontal-), and ENI (Extraterrestrial Irradiance) [W/m²]; albedo [dimensionless],
%       and optionally AMr (relative air mass), and custom Perez coefficients F1, F2.
%
%   SP - structure with fields Az and El (elevation), both Nt-vectors in degrees. Azimuth 
%       convention is assumed to be N2E (astronomical), for SKYREGIONS defined in a system with
%       x = East, y = North.
%
%   model - optional string {'perez','haydavies','igawa'}, determines whether to use the Perez et 
%       al. 1990, Hay & Davies 1980, or Igawa et al. 2004  model. 
%       The default is SimOptions.diffusemodel.
%       To specify a non-standard set of coefficients for the perez model, use the syntax 
%       'perez - SETNAME' where SETNAME can be any of the sub-models in PVLMOD_PEREZ.
%       Alternatively, if fields F1 and F2 are passed along with MD, they will be used as custom
%       anisotropy indices for circumsolar-fraction and horizon-brightening.
%
% OUTPUT:
%   DSKY, DGND, DCS - [Nt,Nc] arrays of radiance values [W/m²sr] for sky-, albedo-, and 
%       circumsolar-regions. 
%
%       (*) NOTE: For Perez and Hay -Davies models, the values are not normalized to the solid-angle
%       of each region, but to the area of their orthogonal projection.
%
% REFERENCES:
% [1] N. Igawa, Y. Koga, T. Matsuzawa, and H. Nakamura, "Models of sky radiance distribution and
%   sky luminance distribution", Solar Energy, vol. 77, no. 2, pp. 137–157, 2004.
% [2] Perez, R., Ineichen, P., Seals, R., Michalsky, J., Stewart, R., 1990. "Modeling daylight 
%   availability and irradiance components from direct and global irradiance. Solar Energy 44 (5).
% [3] Hay, J.E., Davies, J.A., 1980. Calculations of the solar radiation incident on an inclined
%   surface. First Canadian Solar Radiation Data Workshop. Ministry of Supply and Services, Canada.
%
% See also: PVLMOD_PEREZ, PVLMOD_HAYDAVIES, RADIANCE_IGAWA, SHADINGREGIONS, POAIRRADIANCE

    narginchk(4,4);
    if nargin < 4, model = getSimOption('diffusemodel'); end
    if isempty(SR), SR = ShadingRegions(model); end
    
    assert(isa(SR,'ShadingRegions'),'Expecting SHADINGREGIONS object');
    parsestruct(MD,{'BNI','DHI','GHI','ENI','albedo'},'opt',{'AMr','F1','F2'},'-n','-r','-v','-e');
    parsestruct(SP,{'El','Az'},'-n','-f','-v','size',[numel(MD.GHI),1]);

    model = strsplit(lower(model),':'); % take sub-model, e.g. 'perez:france1988'
    switch model{1}
    case 'igawa'
        [Dsky,Dcs] = radiance_igawa(MD,SP,SR);
        
        % PROVISIONAL: isotropic, lambertian albedo
        Dgnd = repmat(MD.GHI.*MD.albedo/pi,1,SR.n.albedo);
        
    case 'perez'
        assert(SR.n.sky == 2 && SR.n.solar == 1 && SR.n.albedo == 1 && SR.haszenith.sky,...
            'SHADINGREGIONS object does not fit the Perez et al. framework');

        if isfield(MD,'F1')
            F1 = MD.F1;
            warning('diffusecomponents:F1','Using custom anisotropy indices');
            if isfield(MD,'F2')
                F2 = MD.F2;
            else
                warning('diffusecomponents:F2','Assuming zero-horizon-brightening, i.e. custom Hay-Davies model');
                F2 = zeros(size(MD.GHI));
            end
        else
            [F1,F2] = pvlmod_perezcoeffs(MD.BNI,MD.DHI,MD.ENI,SP.El,model{2:end});
        end

        prj = polyprojector('ortho','normalize',1);
        [F_skyH,~,F_csH] = viewfactors(SR,prj,[0;0;1],sph2cartV(90-SP.Az,SP.El)');
        [F_skyV,F_gndV] = viewfactors(SR,prj,[1;0;0],[0;0;1]);
        
        if F_skyH(1) < F_skyV(1) || F_skyH(2) > F_skyV(2)
            error('Unexpected sky-region behaviour');
        end
        
        F_skyH = F_skyH.*SR.solidangles.sky;    % projected area of ISO & HB regions
        F_skyV = F_skyV.*SR.solidangles.sky;    % id. (on vertical surface)
        F_gndV = F_gndV.*SR.solidangles.albedo; % projected area of GND on vertical surface

        Dsky(:,2) = F2./(F_skyV(2) - 0.5*F_skyH(2));                           % HB/DHI  (*)
        Dsky(:,1) = (1 - F1 - Dsky(:,2).*F_skyH(2))./(F_skyH(1) + F_skyH(2));  % ISO/DHI (§)
        
        % (*) NOTE: negative F2 values are allowed to effectively increase isotropic component,
        %   by (1-F1-F2·K), but actual horizon band radiance should be clipped to realistic
        %   gradation bounds (CIE I and IV) std skies, with gradation [0.33,7.2]
        GRADBOUNDS = [0.3349,7.178];
        Dsky(:,2) = max(Dsky(:,2), (1-F1)./(F_skyH(1)/GRADBOUNDS(1) + F_skyH(2))-Dsky(:,1));
        Dsky(:,2) = min(Dsky(:,2), (1-F1)./(F_skyH(1)/GRADBOUNDS(2) + F_skyH(2))-Dsky(:,1));
        
        Dsky(:,2) = Dsky(:,2) + Dsky(:,1); % (§) isotropic overlaps HB
        Dsky = Dsky.*MD.DHI;
        
        % Old:
        % Dsky(:,1) = MD.DHI.*(1-F1-F2.*F_skyH(2)/F_skyV(2))./F_skyH(1);
        % Dsky(:,2) = MD.DHI.*max(0.33,F2)/F_skyV(2);

        Acs = SR.solidangles.solar;     
        Dcs = MD.DHI.*F1./(F_csH.*Acs);    % circumsolar
        
        % Old:
        % Acs = pi*sind(SR.cs)^2; % area of projected CS region (AOI = 0)
        % Dcs = MD.DHI.*F1./max(0.087,sind(SP.El))/Acs;
        
        Dgnd = MD.GHI.*MD.albedo./(2*F_gndV);                             % albedo

    case 'haydavies'
        if isfield(MD,'F1')
            F1 = MD.F1;
            warning('diffusecomponents:F1','Using custom anisotropy indices');
        else
            F1 = MD.BNI./MD.ENI; % circumsolar fraction F1 = "anisotropy index" = kn
        end
        if SR.n.sky < 2
            warning('diffusecomponents:F1','Ignoring horizon-brightening indices F2, not enough sky regions');
        end
        
        Dsky = MD.DHI.*(1-F1)/pi;
        Dgnd = MD.GHI.*MD.albedo/pi; % For an isotropic half dome: DHI = pi·radiance
        
        Acs = pi*sind(SR.cs)^2; % area of projected CS region (AOI = 0)
        Dcs = MD.DHI .* F1./max(sind(SP.El),0.01745)/Acs;
    otherwise
        error('Unknown transposition model: %s',strjoin(model,':'));
    end
end