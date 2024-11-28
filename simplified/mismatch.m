function [Fe,opt] = mismatch(Nsb,Ntb,rd,varargin)
% FE = MISMATCH(NSB,NTB,RD,[MDL,X,FF,VD_VMP]) - Estimate mismatch losses according to either 
%   Deline et al. (2013) or Martínez-Moreno et al. (2010).
%
% INPUT:
%   NTB - number of diode-protected cell-blocks (submodules) per string.
%   NSB - number of shaded cell-blocks per string.
%   RD - (1 - shade opacity) = (shaded) diffuse tilted / (non-shaded) global tilted irradiance
%   MDL - 'deline' (default) or 'martinezmoreno'
%   X - fraction of strings that are shaded by NSB (for MDL = 'deline', default = 1)
%   FF - module fill factor (only for MDL = 'deline', default = 0.75)
%   VD_VMP - Nominal diode voltage divided by module MPP voltage: Vd/Vmp, default = 0.01.
%
% OUTPUT: 
%   FE - String output fraction (1 = not shaded, rd = bypassed string)
%
% REF:
% [1] Martínez-Moreno, F., Muñoz, J., Lorenzo, E., 2010. Experimental model to estimate shading 
%   losses on PV arrays. Solar Energy Materials and Solar Cells 94, 2298–2303. 
%   https://doi.org/10.1016/j.solmat.2010.07.029
%
% [2] Deline, C., Dobos, A., Janzou, S., Meydbray, J., Donovan, M., 2013. A simplified model of
%   uniform shading in large photovoltaic arrays. Solar Energy 96, 274–282.
%   https://doi.org/10.1016/j.solener.2013.07.008
%
% See also: INFINITE_ROWS, PV_SIMPLE

    if nargin == 0, test(); return; end

    opt.model = 'deline';
    opt.X = 1;
    opt.FF = 0.75;
    opt.Vd_Vmp = 0.01;
    opt = getpairedoptions(varargin,opt,'dealrest');
    
    opt.model = parselist(opt.model,{'deline','martinezmoreno','worstcase'}');
    
    compatiblesize(Nsb,Ntb,rd,opt.X,opt.FF,opt.Vd_Vmp);
    validateattributes(rd,'numeric',{'real','finite','>=',0,'<=',1});
    validateattributes(opt.X,'numeric',{'real','finite','>=',0,'<=',1});
    validateattributes(opt.FF,'numeric',{'real','finite','>',0.6,'<',0.9});
    validateattributes(Nsb,'numeric',{'real','integer','nonnegative'});
    validateattributes(Ntb,'numeric',{'real','integer','positive'});
    assert(all(Nsb <= Ntb,'all'),'Inconsistent Nsb > Ntb');
    
    switch opt.model
    case 'worstcase'
        Fe = (Nsb == 0).*(1-rd) + rd;
       
    case 'martinezmoreno'
        Fe = martinezmoreno(Nsb,Ntb,rd);
        
    case 'deline'
        S = Nsb./Ntb;
        Fe = deline(opt.X,S,rd,opt.FF,opt.Vd_Vmp);
    end
end

function Fe = martinezmoreno(Nsb,Ntb,rd)
% Approx. Electrical shading losses according to Martínez-Moreno et al. (2010)

    Fe = 1-Nsb./(Ntb+1);    % string missmatch factor (1-Fes)/(1-Fgs)
    Fe = Fe.*(1-rd) + rd;   % apply only to direct fraction
end

function Ps_Ps0 = deline(X,S,Ee,FF,Vd_Vmp)
% Approx. Electrical shading losses according to Deline et al. (2013)
% X = fraction of shaded parallel strings
% S = fraction of shaded submodules for each shaded string
% Ee = shade opacity
% FF = module fill factor

    if nargin < 5 || isempty(Vd_Vmp), Vd_Vmp = 0; end

    C1 = (109*FF - 54.3).*exp(-4.5*X); % eq. 7
    C2 = -6*X.^2 + 5*X + 0.28;
    C2(X > 0.65) = 1;
    Ps1_Ps0 = 1 - C1.*S.^2 - C2.*S;

    Ps2_Ps0 = 1 - S./X.*(1 + Vd_Vmp);

    C3 = (-0.05*Ee-0.01).*X + (0.85*FF - 0.7).*Ee - 0.085*FF + 0.05;
    C3 = max(C3,Ee-1);
    Ps3_Ps0 = C3.*(S -1) + Ee;

    Ps_Ps0 = max(max(Ps1_Ps0,Ps2_Ps0),Ps3_Ps0);
    
    % cases = (Ps_Ps0 == Ps1_Ps0)*1 + (Ps_Ps0 == Ps2_Ps0)*2 + (Ps_Ps0 == Ps3_Ps0)*3; 
    % Psys_Psys0 = X.*Ps_Ps0 + (1-X);
end

function test()

    Ntb = 60;
    Nsb = 0:Ntb;
    rd = (0.2:0.2:0.8)';
    X = [1,1/4];
        
    n = numel(rd);
    m = numel(X);
    
    GUIfigure('mismatch_test','',num2str(m,'%d:1')); clf();
    C = get(gcf,'DefaultAxesColorOrder');
    set(gcf,'DefaultAxesColorOrder',repmat(C(1:n,:),2,1));
    
    for j = 1:m
        subplot(1,m,j); hold on;
    
        Fe = mismatch(Nsb,Ntb,rd,'deline',X(j));
        H{1} = plot(Nsb./Ntb,Fe);

        Fe = mismatch(Nsb,Ntb,rd,'martinezmoreno');
        H{2} = plot(Nsb./Ntb,Fe,'--');

        xlabel('Shading fraction'); 
        ylabel('String output');
        title(sprintf('X = %0.0f%%',X(j)*100));

        if j == 1
            tags = arrayfun(@(e) sprintf('r_d = %0.0f%%',e),rd*100,'unif',0);
            tags(end+1:end+2) = {'Deline et al. (2013)','Martínez-M. et al. (2010)'};

            h = cat(1,H{:},copy(H{1}(1)),copy(H{2}(1)));
            h(end).Color(:) = 0;
            h(end-1).Color(:) = 0;
            legend(h([1:n,end-1:end]),tags,'location','southwest','box','off');
        end
    end
end