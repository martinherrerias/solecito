function [iam,vf] = circumsolarIAM(fIAM,IA,R)
% [IAM,VF] = CIRCUMSOLARIAM(FIAM,IA,CSR) - estimate the integrated IAM for a uniform circumsolar
%   region of radius CSR located at an angle IA from the surface normal.
%
% See also: DIFFUSEIAM, CHECKIAM

    DZ = 1.0;

    if nargin == 0, test(); return; end
    if nargin < 3, R = 25; end
    if nargin < 2, IA = []; end

    if isempty(IA) || numel(IA) > 2*(90+R)/DZ
        ia = unique([linspace(0,90+R,round((90+R)/DZ)+1),180],'sorted'); 
        interp = true;
    else
        ia = min(abs(solarposition.fixtoplusminus180(IA)),90+R);
        interp = false;
    end
    
    fIAM = checkIAM(fIAM);

    TOL = 1e-4;
    warning_reseter = naptime();  %#ok<NASGU>
    warning('off','MATLAB:integral2:maxFunEvalsPass');
    warning('error','MATLAB:integral2:maxFunEvalsFail'); %#ok<CTPCT>
    
    r = @(x) acos(cosd(R)./cosd(x));
    
    b = repmat(R,size(ia));
    a = ia - min(ia + R,90);
    
    msg = {};
    % k = 1./(pi/4*(sin(b)-sin(a)));
    y = zeros(size(ia));
    n = y;
    for j = 1:numel(y)
        if a(j) >= b(j), y(j) = 0; continue; end
        try
            y(j) = integral2(@(x,z) fIAM(ia(j)-x).*cosd(ia(j)-x).*cosd(z),a(j),b(j),0,r,'RelTol',TOL);
            n(j) = integral2(@(x,z) cosd(ia(j)-x).*cosd(z),a(j),b(j),0,r,'RelTol',TOL);
        catch ERR
            msg{end+1} = getReport(ERR); %#ok<AGROW>
        end
    end
    if ~isempty(msg)
       msg = uniquecell(msg);
       warning(strjoin(msg,newline()));
    end
    
    y = y./n;
    n = n/180/(1-cosd(R));
    
    y(ia == 180) = 0;
    
    if ~interp
        iam = y;
        vf = n;
    elseif isempty(IA)
        clip = @(x) abs(solarposition.fixtoplusminus180(x));
        iam = @(x) interp1(ia,y,clip(x));
        vf = @(x) interp1(ia,n,clip(x));
    else
        iam = interp1(ia,y,IA);
        vf = interp1(ia,n,IA);
    end
end

function test()

    R = 30;
    
    [fIAM,iam_txt] = checkIAM('martinruiz','maxloss',0.1);

    x = linspace(50,90+R,101);
    
    c = circumsolarshading(90-x,R);
    c = c.*max(0,cosd(x));
        
    [y,n] = circumsolarIAM(fIAM,x,R);
    b = max(0,cosd(x));
    
    f = fIAM(x);
    
    GUIfigure('circumsolarIAM','circumsolarIAM test'); clf(); hold on;
    
    C = colorcitos(3);

    plot(x,n,'k-','DisplayName','IAM = 1');
    plot(x,y.*n,'k--','DisplayName',strrep(iam_txt,' IAM model',''));
    
    plot(x,n,'color',C(1,:),'DisplayName','Analytic');
    plot(x,b,'color',C(2,:),'DisplayName','Point (center)');
    plot(x,c,'color',C(3,:),'DisplayName',sprintf('Point w/ shading [Varga & Mayer]'));
    legend('box','off','location','northeast');

    plot(x,y.*n,'--','color',C(1,:),'HandleVisibility','off');
    plot(x,cosd(x).*f,'--','color',C(2,:),'HandleVisibility','off');
    plot(x,c.*f,'--','color',C(3,:),'HandleVisibility','off');
    
    plot([1,1]*(90 - R),ylim(),'color',[1 0.5 1],'HandleVisibility','off');
    text(90 - R,0.05, sprintf('  CSR = %0.0f°',R),'color',[1 0.5 1]);
        
    xlabel('Incidence angle');
    ylabel('(View Factor) x IAM');
    grid on;
end 

function c = circumsolarshading(h,r)
% c = circumsolarshading(h,r) - returns the visible fraction of a circle of radius r, h degrees
%   above the horizon.

    f = abs(h) < r;

    b = sqrt(r.^2 - h(f).^2);
    x = 2*atan2(b,h(f));
    
    c = (h > 0)*1;
    c(f) = 1 - (x - sin(x))/(2*pi);
end
