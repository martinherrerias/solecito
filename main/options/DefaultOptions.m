function Defaults = DefaultOptions()
% Returns a structure with default simulation options

    % Inherit defaults from matlab-pvmeteo
    PATH = '../../matlab-pvmeteo';
    Defaults = DefaultsFrom(fileparts(mfilename('fullpath')),PATH);

    Defaults.prjname = '';
    Defaults.version = 'GUI'; % (see runningfromUI)
    Defaults.verbose = true;    
    
%% PROVISIONAL: important user input

    % Auto-completion of missing mount-parameters
    Defaults.mounts.tracklimits.ver = [-150,150,0,90];
    Defaults.mounts.tracklimits.hor = [-45,45];
    Defaults.mounts.axisoffset = [0;0;0];
    Defaults.mounts.backtracking = false;
    Defaults.mounts.groundclearance = 0.5;
    
    % Cell temperature model: pick one of three options:
    Defaults.CellTemp.model = 'SAPM';    % Sandia Cell Temp. model
                                         % Ommit parameters for UI mounting-type pick
                                         
    % Defaults.CellTemp.a_wind = -3.47;  % Sandia Cell Temp. model coefficient 'a'
    % Defaults.CellTemp.b_wind = -.0594; % Sandia Cell Temp. model coefficient 'b'
    % Defaults.CellTemp.delT = 3.0;      % Sandia Cell Temp. model coefficient 'delta T' (K)
    %
    % Defaults.CellTemp.model = 'Uval'; % Faiman quasi-physical model
    % Defaults.CellTemp.Uconst = 29.0;  % Heat transfer coefficient (wind independent) [W/m²K]
    % Defaults.CellTemp.Uwind = 0.0;    % Heat transfer coefficient (wind dependent) [W·s/m³K]
    % Defaults.CellTemp.absort = 0.9;   % optical absorptivity
    % 
    % Defaults.CellTemp.model = 'NOCT'; % Faiman model with Uconst adjusted to NOCT
    % Defaults.CellTemp.NOCT = 44;          % Cell Temperature at NOC (°C)
    % Defaults.CellTemp.Uwind = 0.0;   
    % Defaults.CellTemp.absort = 0.9;
    
    % Perform partial-simulations:
    Defaults.analysedmppts = 'all';     % {'all','trck', or numeric vector}
    
    % Treatment of meteo-data    
    Defaults.meteo.filter = 'wholedays'; % {'all','notdark','wholedays','12days', or MATLAB logical
                                         % expresion, using any field in MeteoData, Time, or SunPos
                                         % structures, plus notdark and wholedays}
    Defaults.meteo.downsample = 1;       % Use e.g. 10 to average 1-min data into 10-min data       
    
%% Shading Analysis / Irradiance Transposition
    Defaults.shadingresolution = 1.25;  % Shading mesh lattice (degrees), 1.25° ~ 5 min
                                        % should be consistent with angtol and meteo-data resolution
	Defaults.diffuseshading = true;     % Otherwise only beam shading analysis
	Defaults.anisotropic = true;        % Include Circumsolar & Horizon Brightening (recommended)
	Defaults.groundshading = false;     % Consider albedo shading by tracker shades on the ground (slow!)
    Defaults.groundshadeangle = 6.0;    % Stop calculating ground shades for sunel < x
	Defaults.flathorizon = false;       % Ignore horizon profile(s)

    Defaults.diffusemodel = 'igawa';    % {'perez','haydavies','igawa'}
    Defaults.skyregions = 'dumortier';  % Any valid input to SHADINGREGIONS e.g. {'dumortier','perez','sunto','eko','bucky','tinN','rbN','urN',...}
    Defaults.HBBwidth = 6.5;            % Perez model Horizon Brightening width (degrees)
	Defaults.CSradius = 25;             % Perez model Circumsolar Radius (degrees)
    Defaults.analysedpts = [-0.4,-1.1]; % [n,m] Use QUADRATUREPOINTS([m,n],A,B) as diffuse-shading 
                                        % analysis points, where A,B are the mount's rectangular 
                                        % limits (A = -[W,H]/2, B = [W,H]/2), and ...
                                        %   if m,n > 0, use m rows, n columns of points
                                        %   if m,n < 0, use ceil(H·m), rows and ceil(W·n) columns 
	
    Defaults.diode = 'schottky';   % 'pn' (~0.6 V), 'schottky' (~0.4 V) or 'smart' (~50 mV), see BYPASSDIODEFIT
    Defaults.material = 'multisi'; % used by PVL_FSSPECCORR for spectral correction, unless ODM
                                   % includes recognized 'material' field
    
%% Electrical Solver
    Defaults.bypassdiodes = true;       % Bypass diodes in parallel with every cell-block inside the module
	Defaults.modulediodes = false;      % An additional bypass diodes in parallel with the complete module (unusual)
	Defaults.blockingdiodes = false;    % Diodes in series with every string (unusual)
	Defaults.parallelelements = false;  % Are cell blocks in parallel inside every module?
    Defaults.autosavemin = 20;          % Save workspace every x minutes
    
    Defaults.celldefects = false;       % requires manual modeling of cell-defects
    
    % Output
    Defaults.solverlog = false;
    Defaults.resultsxls = false;
    
%% Terrain
    Defaults.maxlocdist = 5;               % (km) Issue warning when meteo data is more than maxlocist away from project
    Defaults.terrain.demradius = 25000;    % (meters) Recommended minimum span of DEM around site (for far-horizon)
    Defaults.terrain.tilethreshold = 0.05; % Don't demand DEM tiles that represent less than this fraction of a demradius-circle
    Defaults.terrain.minringwidth = 10;    % (DEM-grid-lattices) Used in GETHORIZONS to set the 'near-ground' limit
    Defaults.terrain.buffer = 50;          % (meters) Used in READTERRAIN to smooth-out DEM errors
    Defaults.terrain.flatthreshold = 0.5;  % (Degrees) Consider a horizon with max(el) < flatthreshold as completely flat
    Defaults.terrain.minmeshdist = 0.5;    % (meters) min. terrain mesh distance (also used for ground shading)
    Defaults.terrain.maxmeshdist = 10.0;   % (meters) max.terrain mesh edge distance
    Defaults.terrain.meshtol = 0.2;        % (meters) allowable terrain-mesh simplification error
    
    Defaults.terrain.offsetcheck = true;
    Defaults.terrain.offsetradius = @(r) max(100,2*r);
    Defaults.terrain.offsetlattice = 40;
    Defaults.terrain.offsetwarning = @(g) max(1,0.1*g);
    Defaults.terrain.offsetP = 0.95;
                                        
%% Solver
    Defaults.MPPmethod = 'IVpp';        % 'interp','ivpp','odm' for uniform cases
	Defaults.savingIVpps = false;       % save all array IVpp's (memory!)
    Defaults.quickshading = true;       % Cell shading factors based on the number of shaded cell-border vertices
    Defaults.skipnonshaded = true;      % Simplify calculation for not-shaded time-steps (slow!)
                                        % consider reducing Gvartol before setting this to false.
%% Precision and Tolerances
    Defaults.cellshadingthreshold = 1/10; 
    
    % angtol: (Degrees) Controls polygon detail (for representation of curves on the unit sphere)
    %   and angular tolerances in general. Recommended is 360/ceil(2·pi/sqrt(6·RelTol)):
    % 	The relaive error in the area of a circle approximated by n segments is 
    %       e = 1-sin(s)/s = s^2/3! - s^4/5! + ... ~ s²/6, for s = 2·pi/n << 1
    %   Taking the next-lower integer fraction of 360° allows expressions like 0:angtol:360, e.g.
    %   360/256 ~ sqrt(6·1e-4)·180/pi for a RelTol of 1e-4.
	Defaults.angtol = 360/ceil(2*pi/sqrt(6*Defaults.RelTol));

    % masktolerance: Maximum fraction of visible field occupied by an obstacle that can be ignored
    % during the beam-shading calculation, pi·(16/60·pi/180)² / 2·pi ~ 1e-5 is the solar disc.
    % For diffuse-shading, masks are calculated from RelTol.
	Defaults.masktolerance = 1e-5;       
    
    % cellshadingtol: If ~quickshading, cell-shading factors are rounded to multiples of this value
    Defaults.cellshadingtol = Defaults.RelTol;     
    
    % Gvartol: Irradiance variability within a string, above which circuit must be solved
    % Clustering (dg/g)² < RelTol ensures that Ns·Pmp(Gavg)/sum(Pmp(G)) - 1 ~ RelTol (#)
    Defaults.Gvartol = sqrt(Defaults.RelTol);            
    
    % (#) NOTE: there is an error in MPP when slightly changing effective irradiance G, aside from 
    % the linear current-scaling Pmp(k·G) ~ k·Pmp(G). This is due to a shift in true MPP-voltage.
    % Theoretically, for k ~ 1.0, 1 - P(Vmp0,k·G)/k·P(Vmp0,G) = (1/4)·(1-k)²
    % That is, the relative error in Pmp that arises from assuming a linear Pmp(k·G) ~ k·Pmp(G) is
    % bounded by 0.25·(1-k)², where |1-k| << 1 is the relative change in effective irradiance.
    % For a given Pmp relative tolerance t, this implies that |1-k| < sqrt(4·t)
    % PENDING: When combining several mismatched modules in series, calculations suggest that the  
    % error is closer to sqrt(dG/G)

%% ODM Interpolant Construction
	Defaults.Lims_Voc0 = [-50,1.5];     % Limits for IV curve (voltage in Voc0 units)
    Defaults.Lims_Isc0 = [-2,2];        % Limits for IV curve (voltage in Isc0 units)
	Defaults.Gegrid = 0:200:1600;       % Grid of irradiances for evaluation of ODM, recommended gap ~ 200 W/m² for RelTol= 1e-4 (will be refined)
	Defaults.Tcgrid = -32:8:96;         % Grid of cell temperatures [...], recommended gap ~ 4°C for RelTol= 1e-4 (will be refined)
	Defaults.Tagrid = -16:4:48;         % Grid of ambient temperatures for Cell-Temp. interpolant
	Defaults.vwgrid = [0,16];           % Grid of wind speeds [...]
    Defaults.meshref = 5;               % Refining factor for scalar (Pmp,Imp,Isc...) interpolants
    
% PROVISIONAL: pre-processor issues
    Defaults.fliptilt = false;          % change if tilt convention doesn't follow PP-specs 2.0
end
