% Load a diode- and module-interpolants
basedir = cd(fileparts(mfilename('fullpath')));
cd('../../../../tests/ODM');

[~,M] = ODM_testbase();
D = BypassInterpolant();
M = M(randi(numel(M))); % keep one of the modules at random

cd(basedir);

% Calculate cell-temperature function
celltemp = getcelltempfcn(M.source,@M.MPPefficiency);

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

fprintf('Input overhead...'); tic()

block_sizeset = cat(1, block_curves.n)';

block_pwl = [cat(1, block_curves.x) cat(1, block_curves.y)];

diode_sizeset = cat(1, diode_curves.n)';

diode_pwl = [cat(1, diode_curves.x) cat(1, diode_curves.y)];

%{
block_pwl = [];
block_sizeset = [];

for j = 1:Ns
    for i = 1:Nm
%        disp([block_curves(i,j).x,block_curves(i,j).y]);
        block_pwl = [block_pwl;[block_curves(i,j).x,block_curves(i,j).y]];
        block_sizeset = [block_sizeset,block_curves(i,j).n];
    end
end

diode_pwl = [];
diode_sizeset = [];
for j = 1:Ns
    for i = 1:Nm
%        disp([block_curves(i,j).x,block_curves(i,j).y]);
        diode_pwl = [diode_pwl;[diode_curves(i,j).x,diode_curves(i,j).y]];
        diode_sizeset = [diode_sizeset,diode_curves(i,j).n];
    end
end
%}
toc()

%disp(block_pwl);
fprintf('The benchmark starts... '); 
fprintf('\n\n');
[out1, out2]=mexcppmdpwltest(block_pwl, diode_pwl, block_sizeset, diode_sizeset);

line(out1', out2');
hold on;
vv = linspace(0,Nm*M.Voc0,100);
axis([-2*Nm Nm*M.Voc0 0 Ns*M.Isc0])
xlabel('Array Voltage [V]');
ylabel('Array Current [A]');
grid on




