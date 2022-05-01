
% Load a diode- and module-interpolants
M = load('test_interpolants.mat','ModIVint','Diode');
D = M.Diode;
M = M.ModIVint;

% Calculate cell-temperature function
celltemp = getcelltempfcn(M.source,'eff',@M.MPPefficiency);
%%

Nm = 20;  % Modules per string
Ns = 6;   % Strings per MPPT

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

fprintf('Adding bypass diodes... '); tic()
module_curves = arrayfun(@(m,d) addparallel([m,d],[-Inf Inf -M.Isc0 2*M.Isc0]),block_curves,diode_curves);
toc()

fprintf('Getting isolated MPPs...'); tic()
Pmp0 = arrayfun(@mpp,module_curves);
toc()

% Add module curves in series to form strings
fprintf('Adding blocks in series...'); tic()
string_curves = arrayfun(@(j) addseries(module_curves(:,j)),1:Ns);
toc()

% Add strings in parallel to get array curve
fprintf('Adding strings in parallel...'); tic()
array_curve = addparallel(string_curves,[-Inf Inf 0 M.Isc0*Nm]);
toc()

fprintf('Getting array MPP...'); tic()
[~,Vmp] = array_curve.mpp();
Imp = array_curve.val(Vmp);
toc()

% Plot
clf()

subplot(1,3,2);
c = jet(Ns);
set(gca,'colororder',c);
plot(string_curves);
axis([-2*Nm Nm*M.Voc0 0 M.Isc0])
xlabel('String Voltage [V]');
ylabel('String Current [A]');
grid on

subplot(1,3,1);
c = c + (rand(Ns,3,Nm)*2-1)*0.2;
c = min(max(0,c),1);
c = permute(c,[2 3 1]);
c = reshape(c,3,[])';
set(gca,'colororder',c);
plot(module_curves);
axis([-2 M.Voc0 0 M.Isc0])
xlabel('Module Voltage [V]');
ylabel('Module Current [A]');
grid on

subplot(1,3,3);
plot(array_curve)
hold on;
plot(Vmp,Imp,'go');
vv = linspace(0,Nm*M.Voc0,100);
plot(vv,sum(Pmp0(:))./vv,':');
plot(scale(M.getIVpp(mean(G(:)),mean(Tc(:))),Nm,Ns))
axis([-2*Nm Nm*M.Voc0 0 Ns*M.Isc0])
xlabel('Array Voltage [V]');
ylabel('Array Current [A]');
grid on
legend('result','MPP','Pmp0','avg')
