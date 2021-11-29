% Load a diode- and module-interpolants
M = load('test_interpolants.mat','ModIVint','Diode');
D = M.Diode;
M = M.ModIVint;

% Calculate cell-temperature function
celltemp = getcelltempfcn(M.source,@M.MPPefficiency);
%%

Nm = 20;  % Modules per string 20
Ns = 6;   % Strings per MPPT 6

Gm = 800; % Mean global effective irradiance
Dm = 200; % Mean diffuse irradiance

Ta = 25;  % Ambient temperature
ws = 2.0; % wind speed

Ss = 50;  % irradiance std. among strings
Sm = 5;   % irradiance std. among modules of the same string
Ps = 0.3; % probability that one string is shaded
Pm = 0.8; % probability that one module is shaded, given that the string is shaded

G = randn([Nm,Ns])*Sm+randn(1,Ns)*Ss + Gm - (rand([1,Ns]) < Ps).*(rand([Nm,Ns]) < Pm)*(Gm-Dm);
Tc = celltemp(G,repmat(Ta,[Nm,Ns]),repmat(ws,[Nm,Ns]));

fprintf('\n\n');

% Get module & bypass-diode curves, add them in parallel
fprintf('Interpolating (%d) curves... ',numel(G)); tic()
block_curves = M.getIVpp(G,Tc);
diode_curves = scale(D.getIVpp(Tc),3,1);
toc()

%disp(block_pwl);
fprintf('The benchmark starts... '); 
fprintf('\n\n');
mexcppmdpwltest_obj(block_curves, diode_curves);





