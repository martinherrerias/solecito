function mexaddseriesparallel()

% Load a diode- and module-interpolants
M = load('test_interpolants.mat','ModIVint','Diode');
D = M.Diode;
M = M.ModIVint;

% Calculate cell-temperature function
celltemp = getcelltempfcn(M.source,'eff',@M.MPPefficiency);
%

for k = 1:10

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

% copy the resource file to the current path
psrc = '../pvProj/Resources';
pdest = '.';
filename = 'SimOption.xml';
source = fullfile(psrc,filename);
dest = fullfile(pdest,filename);
copyfile(source, dest);

tic()
for j = Ns:-1:1
    XY = arrayfun(@(p) [p.x,p.y],block_curves(:,j),'unif',0);
    [x, y, size] = mexaddseries(Nm,XY{:},[0,0,0],[-Inf Inf -Inf Inf]);
    string_curves(j) = mdpwl(x(1:size)', y(1:size)',0);
end
XY = arrayfun(@(p) [p.x,p.y],string_curves,'unif',0);
[x, y, size] = mexaddparallel(Ns,XY{:},[0,0,0],[-Inf Inf -Inf Inf]);
array_curve = mdpwl(x(1:size)', y(1:size)',0);
toc()

tic()
for j = Ns:-1:1
    string_curves(j) = addseries(block_curves(:,j),[-Inf Inf -Inf Inf],[0,0,0]);
end
array_curve = addparallel(string_curves,[-Inf Inf -Inf Inf],[0,0,0]);
toc()
   
end


%%
% 
% fprintf('Input overhead...'); tic()
% 
% block_sizeset = cat(1, block_curves.n)';
% 
% block_pwl = [cat(1, block_curves.x) cat(1, block_curves.y)];
% 
% diode_sizeset = cat(1, diode_curves.n)';
% 
% diode_pwl = [cat(1, diode_curves.x) cat(1, diode_curves.y)];
% 
% %{
% block_pwl = [];
% block_sizeset = [];
% 
% for j = 1:Ns
%     for i = 1:Nm
% %        disp([block_curves(i,j).x,block_curves(i,j).y]);
%         block_pwl = [block_pwl;[block_curves(i,j).x,block_curves(i,j).y]];
%         block_sizeset = [block_sizeset,block_curves(i,j).n];
%     end
% end
% 
% diode_pwl = [];
% diode_sizeset = [];
% for j = 1:Ns
%     for i = 1:Nm
% %        disp([block_curves(i,j).x,block_curves(i,j).y]);
%         diode_pwl = [diode_pwl;[diode_curves(i,j).x,diode_curves(i,j).y]];
%         diode_sizeset = [diode_sizeset,diode_curves(i,j).n];
%     end
% end
% %}
% toc()
% 
% %disp(block_pwl);
% fprintf('The benchmark starts... '); 
% fprintf('\n\n');
% [x, y]=mexcppmdpwltest(block_pwl, diode_pwl, block_sizeset, diode_sizeset);
% 
% line(x', y');
% hold on;
% vv = linspace(0,Nm*M.Voc0,100);
% axis([-2*Nm Nm*M.Voc0 0 Ns*M.Isc0])
% xlabel('Array Voltage [V]');
% ylabel('Array Current [A]');
% grid on




