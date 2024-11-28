function x = circumsolarshading(varargin)
% X = CIRCUMSOLARSHADING(H,R) - returns the visible fraction of a circle of radius R,
%   when its center is H (units of R) distant from an opaque semiplane. -R < H < 0 implies the
%   center of the circle is included in the semiplane, and a distance H from its edge.
%
%   For infinite-row shading, H is what Varga & Mayer [1] call (alpha_s' - Psi), the difference
%   between the projected solar elevation and the shade free angle:
%
%       SFA = INFINITE_ROWS(pitch,table_width,tilt,[sections]);  % Psi
%       prjel = atan2d(tand(sunel),cosd(sunaz-surfaz));          % alpha_s'    
%       C = CIRCUMSOLARSHADING(prjel - SFA,R);
%
% X = CIRCUMSOLARSHADING(H1,H2,PHI,R) - returns the visible fraction of a circle of radius R,
%   when its center is a distance H1 from one semiplane, and H2 from another, both semiplanes 
%   being at an angle PHI from each other.
%
% [1] Varga, N., Mayer, M.J., 2021. Model-based analysis of shading losses in ground-mounted 
%   photovoltaic power plants. Solar Energy 216, 428–438. 
%   https://doi.org/10.1016/j.solener.2021.01.047
%
% See also: INFINITE_ROWS

    if nargin == 0, test(); return; end
    narginchk(2,5);

    if ischar(varargin{end}) && strcmpi(varargin{end},'-deg')
        uses_degrees = true;
        varargin(end) = [];
    else
        uses_degrees = false;
    end
    
    for j = 1:numel(varargin)-1
        validateattributes(varargin{j},'numeric',{'real'});
        if uses_degrees
            varargin{j} = solarposition.fixtoplusminus180(varargin{j});
        end
    end
    validateattributes(varargin{end},'numeric',{'real','positive'});
    [varargin{1:end}] = compatiblesize(varargin{1:end});

    switch numel(varargin)
        case 2, x = simpleshading(varargin{:});
        case 4, x = composedshading(varargin{:});
        otherwise, error('Expecting two (H,R) or four (H1,H2,PHI,R) arguments');
    end
end

function x = simpleshading(h,r)
% NOTE: there seems to be an error in [1], their eq. (16) is equivalent to (X - sin(X/2))/(2pi)
% vs the correct (X - sin(X))/(2pi), where X = 2*acos(h/r), as in their diagram.

    f = abs(h) < r; % partial shading

    X = 2*acos(h(f)./r(f));

    x = (h > 0)*1;
    x(f) = 1 - (X - sin(X))/(2*pi);    
end

function [x,varargout] = composedshading(h1,h2,phi,r)

    sz = size(h1);
    x = zeros(sz);
    
    % 1. Fully Shaded by either h1 or h2
    fullyshaded = (h1 <= -r) | (h2 <= -r);
    x(fullyshaded) = 0;
    
    % td = To Do, cc = Current Case
    td = ~fullyshaded;

    % 2a. Only shaded by h2, if at all
    cc = td & (h1 >= r);
    x(cc) = simpleshading(h2(cc),r(cc)); 
    
    td(cc) = false;
    
    % 2b. Only shaded by h1, if at all
    cc = td & (h2 >= r);
    x(cc) = simpleshading(h1(cc),r(cc));
    
    td(cc) = false;
    b = td; % Base index, to calculate A1, A2, alpha, ...

    th1 = acos(h1(b)./r(b));
    th2 = acos(h2(b)./r(b));
    
    phi = mod(abs(phi(b)),360);
    phi = min(phi,360-phi)*pi/180;
    cphi = cos(phi);
    sphi = sin(phi);
    
    c1m2 = cos(th1-th2);
    c1p2 = cos(th1+th2);
    
    alpha = (th1 + th2 - phi)/2;
    y = sqrt(h1(b).^2 + h2(b).^2 - 2*h1(b).*h2(b).*cphi)./sphi;

    A1 = 0.5*(2*th1 - sin(2*th1));
    A2 = 0.5*(2*th2 - sin(2*th2));

    % 3. Separate shadows (no intersection)
    cc = b;
    cc(b) = td(b) & (alpha < 0);
    j = cc(b);
    x(cc) = 1 - (A1(j)+A2(j))/pi;
    
    td(cc) = false;
    
    % 4. One of the shadows fully contains the other
    cc = b;
    cc(b) = td(b) & (alpha > th1) | (alpha > th2);
    j = cc(b);
    x(cc) = 1 - max(A1(j),A2(j))/pi;
    
    td(cc) = false;
    
    % 5. Fully shaded by compound shadow
    cc = b;
    cc(b) = td(b) & (abs(y) > r(b));
    x(cc) = 0;
   
    td(cc) = false;
    
    % 6. Compound shadow
    j = td(b);
    x(td) = 1 - (A1(j) + A2(j) + phi(j) - sphi(j).*c1m2(j) + ...
            (cphi(j) - c1p2(j)).*(1-cphi(j).*c1m2(j))./sphi(j))/(2*pi);
        
    if nargout > 1
        varargout = repmat({NaN(sz)},1,5);
        varargout{1}(b) = th1*180/pi;
        varargout{2}(b) = th2*180/pi;
        varargout{3}(b) = phi*180/pi;
        varargout{4}(b) = alpha*180/pi;
        varargout{5}(b) = y;
    end
end

function test()

    N = 50; % batch size
    M = 10; % batches
    
    c = polygon(360);

    for k = 1:M
        
        h = rand(N,2)*2.2-1.1;
        th = rand(N,2)*180 + rand(N,1)*180;
            
        x0 = zeros(N,1);
        for j = 1:N
            x0(j) = calcwithpolygons(h(j,:),th(j,:));
        end
        
        [x,th1,th2,phi,alpha] = composedshading(-h(:,1),-h(:,2),diff(th,1,2),1);
           
        bad = (x0 - x) > 0.01;
        if any(bad), j = find(bad,1);
        else, j = N;
        end
        [~,p,r] = calcwithpolygons(h(j,:),th(j,:));
        
        GUIfigure('circumsolarshading'); clf(); hold on;
        polyplot(c);
        polyplot(p,'none');
        polyplot(r);
        axis equal; axis([-1,1,-1,1]); set(gca,'visible',false);

        txt = sprintf(['x_0 = %0.1f%%\nx = %0.1f%%\n\n\\theta_1 = %0.1f\n',...
                       '\\theta_2 = %0.1f\n\\phi = %0.1f\n\\alpha = %0.1f'],...
                       x0(j)*100,x(j)*100,th1(j),th2(j),phi(j),alpha(j));
        text(0.3,0.3,txt,'interpreter','tex','fontsize',14,'verticalalignment','top');
        
        if any(bad), keyboard(); end
        pause();
        % keyboard();
    end
    
    function [f,p,r] = calcwithpolygons(H,TH)
        p = arrayfun(@(h,th) polyrotate(polygon([-2,2],[-2 h]),th),H,TH);
        % p = polyrotate(polygon([-1,1],[-1 h(1)]),th(1));

        r = substractpolygons(c,p,'pos','pos');
        f = r.area/c.area;
    end
end
