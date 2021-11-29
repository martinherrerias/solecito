function varargout = diffuseIAM(fIAM,tilt)
% [ISO,ALB,HB] = DIFFUSEIAM(FIAM) - Integrate IAM function under an isotropic radiance & flat
%   terrain assumption, returning function handles for the tilted components:
%
%   ISO(t) - Isotropic component IAM as a function of surface tilt t. Use ISO(t+a) to consider
%       infinite-row shading with shade-free angle a.
%   ALB(t) - Albedo component IAM as a function of surface tilt t
%   HB(t) - Horizon-brightening component IAM as a function of surface tilt t
%
%   TODO: integral2 crashes with discontinuous functions, e.g. checkIAM('hcpv')
%
% See also: CHECKIAM

    DZ = 1.0; % degrees, resolution of interpolation structure for f.iso,...
    
    fIAM = checkIAM(fIAM);

    df = @(x,z) cos(x).*cos(z).*fIAM(acosd(cos(x).*cos(z)));    % cos(ia)·fIAM(ia)
    lineIAM = @(z) integral(@(x) df(x,z),0,pi/2)/cos(z);

    
    if nargin < 2 || numel(tilt) > 2*(90/DZ)
    % ... pre-calculate an interpolation set, for performance
        zz  = linspace(0,pi/2,ceil(90/DZ));
        interp = true;
    else
        zz = (90 - tilt)*pi/180;
        interp = false;
    end
    
    iam_iso = wedgeIAM(-pi/2,zz); 
    iam_alb = wedgeIAM(zz,pi/2); 
    iam_hb = arrayfun(@(z) lineIAM(z),zz);
    
    if ~interp, varargout = {iam_iso,iam_alb,iam_hb}; return; end
    
    if nargin < 2
        iso = @(T) interp1(zz,iam_iso,(90 - T)*pi/180);
        alb = @(T) interp1(zz,iam_alb,(90 - T)*pi/180);
        hb = @(T) interp1(zz,iam_hb,(90 - T)*pi/180);
    else
        iso = interp1(zz,iam_iso,(90 - tilt)*pi/180);
        alb = interp1(zz,iam_alb,(90 - tilt)*pi/180);
        hb = interp1(zz,iam_hb,(90 - tilt)*pi/180);  
    end
    varargout = {iso,alb,hb};
    
    function r = wedgeIAM(a,b)
    % Equivalent to:
    %   r = arrayfun(@(a,b) integral2(@(x,z) cos(x).*df(x,z),0,pi/2,a,b)/(pi/4*(sin(b)-sin(a)));
    %
    % Provisionally solves integration convergence issues for discontinuous functions -- e.g.
    % diffuseIAM('hcpv') -- by lowering tolerance, and capturing warnings.
    %
    % TODO: find robust integration solution for discontinuous functions
    % wedgeIAM = @(a,b) integral2(@(x,z) cos(x).*df(x,z),0,pi/2,a,b)/(pi/4*(sin(b)-sin(a)));

        TOL = 1e-4;
        warning_reseter = naptime();  %#ok<NASGU>
        warning('off','MATLAB:integral2:maxFunEvalsPass');
        warning('error','MATLAB:integral2:maxFunEvalsFail'); %#ok<CTPCT>
        
        [a,b] = compatiblesize(a,b);
        msg = {};
        k = 1./(pi/4*(sin(b)-sin(a)));
        r = NaN(size(k));
        for j = 1:numel(a)
            if a(j) >= b(j), r(j) = 0; continue; end
            try
                r(j) = k(j).*integral2(@(x,z) cos(x).*df(x,z),0,pi/2,a(j),b(j),...
                    'AbsTol',TOL*k(j),'RelTol',TOL);
            catch ERR
                msg{end+1} = getReport(ERR); %#ok<AGROW>
            end
        end
        if ~isempty(msg)
           msg = uniquecell(msg);
           warning(strjoin(msg,newline()));
        end
    end
end

