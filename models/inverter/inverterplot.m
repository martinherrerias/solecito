function varargout = inverterplot(Inverter)
% INVERTERPLOT(Inverter) - Print inverter efficiency curves
% See also: GUIMODELS, ONDREAD, ONDEVAL, PVL_SNLINVERTER

T0 = 25;

p = linspace(Inverter.Ps0,Inverter.Pdc0,1000);
v = linspace(Inverter.MPPTLow,Inverter.MPPTHi,3);

[P,V] = ndgrid(p,v);
etaV = ONDeval(Inverter,P,V,T0*ones(size(P)))./P*100;
legstr = arrayfun(@(v) sprintf('%0.0f V, %d°C',v,T0),v,'unif',0);
%lim = [prctile(etaV(:),20),max(etaV(:))];

if isfield(Inverter,'TPLim') && Inverter.TPLim(0) ~= Inverter.TPLim(100)
    ta = Inverter.TPLim.GridVectors{1}(2:5);
    [~,idx] = unique(Inverter.TPLim(ta),'stable');
    ta = ta(idx);
    [P,~] = ndgrid(p,ta);
    etaT = min(Inverter.TPLim(ta)'./P*100,100);
    ta = ta(any(etaT > 0,1));
    etaT = etaT(:,any(etaT > 0,1));
    legstr = cat(1,legstr',arrayfun(@(t) sprintf('\\eta_{max} @ %0.0f°C',t),ta,'unif',0));
    %lim(1) = min(lim(1),prctile(etaT(:),50));
else, etaT = zeros(0,0);
end

figh = GUIfigure('inverter','Inverter Model'); clf();
% set(figh,'defaultaxescolororder',jet(numel(legstr)));

if nargout > 0, varargout{1} = figh; end

hold all
arrayfun(@(v) plot(p/Inverter.Pdc0,etaV(:,v)),1:size(etaV,2));
arrayfun(@(t) plot(p/Inverter.Pdc0,etaT(:,t)),1:size(etaT,2));
hold off

yy = [80:2:90 91:100];
lbl = arrayfun(@num2str,yy,'unif',0);
lbl(yy > 90 & mod(yy,2)>0) = {''};
% yy = unique(round(1-logspace(log10(1-lim(1)/100),log10(1-lim(2)/100),9),3)*100);
% yy(end+1) = min(100,round(yy(end) + 0.05*(yy(end)-yy(1)),1));
% yy = unique(yy);
yticks(yy);
yticklabels(lbl);
axis([0 1 yy(1) yy(end)]);
% axis square
grid on

legend(legstr(:),'fontsize',8,'location','sw');
legend boxoff
title(Inverter.name); set(get(gca,'title'),'interpreter','none')
xlabel(sprintf('Input Power / %0.2f kW',Inverter.Pdc0/1000)); 
ylabel('Efficiency [%]');

