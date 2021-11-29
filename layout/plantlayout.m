function [F,Vfcn,sunvec,R] = plantlayout(Trck,sunaz,sunel,AzConv)
% [F,V,S,R] = PLANTLAYOUT(TRCK,SUNAZ,SUNEL,[AZCONV]) - Calculate the geometry of a PV plant
%   with layout TRCK, at solar position SUNEL, SUNAZ (in azimuth convention AZCONV).
%
%   TRCK - structure with fields {centers,origin,rotation} as required by CHECKCOORDSYSTEM;
%       {type, azimuth, tilt, slope, tracklimits, backtracking, groundcoverratio} as required
%       by MOUNTROTATIONS;  TRCK.geom (a pvArea object with TRCK.geom.border representing
%       the mount outline), and an optional TRCK.axisoffset (see MOUNTROTATIONS).
%
%   F, V - Index array F [Nu,p] and list of vertices V [Nu·p,3] are a face-vertex representation 
%       of the Nu rotated and translated mounts in TRCK. Vertices V are in the Project Coordinate 
%       System (PCS - see MOUNTROTATIONS for notes on coordinate system conventions)
%
%   S - [Nt,3] (stack of) unit vector(s) in PCS pointing towards the sun.
%   R - [3,3,Nu*,Nt*] stack of rotation matrices (one for each mount and timestep, with 
%       dimensions collapsed for static systems and/or systems with common rotation).
%
%  [F,VFCN,S,R] = PLANTLAYOUT(...) - if SUNAZ, SUNEL are not scalars, and if the mounts are not
%       static, a function handle VFCN will be returned instead of the array of vertices V. The
%       (changing) set of vertices for time t (1..Nt) is then to be retrieved by V = VFCN(t).     
%
%  TODO: this should be part of a class, where F and V are just dependent properties or methods
%       that take SUNEL, SUNAZ.
%
% See also: MOUNTROTATIONS, CHECKCOORDSYSTEM, GROUNDSHADING, PLOTTRACKERARRAY

    if nargin < 2, sunaz = []; end
    if nargin < 3, sunel = []; end
    if nargin < 4 || isempty(AzConv), AzConv = 'N2E'; end
    if ~isequal(AzConv,'N2E')
        sunaz = solarposition.fixazimuth(sunaz,AzConv,'N2E',Trck.origin(2));
    end
    
    parsestruct(Trck,{'centers'},'numeric','real','nonempty','size',[3,NaN]);
    parsestruct(Trck,{'geom'},'class','pvArea','scalar',@(x) isscalar(x.border));
    parsestruct(Trck,{'origin'},'numeric','real','size',[1,3]);
    Ntr = size(Trck.centers,2);

    % Let mountrotations parse the rest of Trck, sunaz and sunel
    R = mountrotations(Trck,sunaz,sunel,'N2E');   % 3·3·Ntr*·Nt* array
        
    if isempty(sunaz) || isempty(sunel) % would have crashed above if not static
        sunvec = NaN(1,3);
    else
        sunaz = 90 - sunaz - Trck.rotation; % project convention: X = 0, Y = 90
        sunvec = sph2cartV(sunaz,sunel);
    end

    % Tracker at starting position
    [v0,~,F] = poly2vef(Trck.geom.border,1);
    if isempty(v0), v0 = [0,0]; end
    
    v0(:,3) = 0;
    if isfield(Trck,'axisoffset')
        validateattributes(Trck.axisoffset,{'numeric'},{'vector','finite','real','numel',3});
        v0 = v0 + Trck.axisoffset(:)';
    end
    % Overall we'll have Ntr faces (trackers) with size(v0,1) vertices each
    F = F + size(v0,1)*(0:Ntr-1)';
    
    if rank(v0) < 3
        
    else
        
    end

    V0 = [];
    if size(R,4) == 1
        V0 = rotatetrackers(1);
        Vfcn = V0;
    else
        Vfcn =  @rotatetrackers;
    end

    function [V,s] = rotatetrackers(t)

        nv = size(v0,1);
        if size(R,4) > 1 || isempty(V0)
        % These are all forms of writing: V = R·v0 + C
            V = repmat(v0,Ntr,1,1);
            if size(R,3) == 1
                V = V*R(:,:,1,t)' + repelem(Trck.centers',nv,1,1);
            else
                V = sum(repelem(permute(R(:,:,:,t),[3,1,2]),nv,1,1).*permute(V,[1,3,2]),3);
                V = V + repelem(Trck.centers',nv,1,1);
            end
        else
            V = V0;
        end
        s = sunvec(t,:)';
    end
end
