function [f,info,IAM,warnmsg] = checkIAM(IAM,varargin)   
% [FIAM,INFO] = CHECKIAM(IAM) - Takes a structure/function-handle IAM, defining an Irradiance-
%   Incidence-Angle modifier model IAM, and returns a (parsed) function handle F, along with an
%   optional INFO structure.
%
%   CHECKIAM() - search the base workspace for a pre-parsed 'IAM' structure or function handle, 
%       that is, one that passes checkIAM(IAM) without warnings. Otherwise use SimOptions.IAM.
%   CHECKIAM(1) - use a perfect absorber, IAM = @(x) 1·(x < 90)
%   CHECKIAM('model') - use specified model, with default parameters
%   CHECKIAM(S) - use (incomplete) parameter structure S
%
% [FIAM,INFO] = CHECKIAM(IAM,'-plot') - plot the IAM function on GUIfigure('iam',...), with
%   the resulting integrated functions F(q) for an isotropic sky & albedo and a horizon-line,
%   for a surface with tilt-angle q, assuming a flat horizon.
%
% INPUT: 
%   If IAM is a structure, depending on IAM.model {'none','ashrae','refraction','hcpv'}, 
%   additional fields are required [IAM.model (required fields) : notes]
%     
%     'ashrae' (maxIncAngle / b) : see PVL_IAM_ASHRAE
%     'refraction' (refractionIndex, glazingThickness, extintionCoeff): see PVL_IAM_PHYSICAL
%     'martinruiz' (ar / maxloss*) : see PVL_IAM_MARTINRUIZ
%     'hcpv' (maxIncAngle) : sharp threshold, f = @(x) 1·(x <= maxIncAngle)
%     'none' () : perfect absorber, f = @(x) 1.0
% 
%   IAM can also be a function handle f : [0,180] → [0,1] that takes the incidence angle(s) in
%     degrees and returns the absorbed fraction(s) of incoming radiation.
%
%   IAM = 'name' is equivalent to IAM = struct('model','name')
%
% OUTPUT:
%   FIAM - function handle FIAM : [0,180] → [0,1] that takes the incidence angle(s) in degrees
%     and returns the absorbed fraction(s) of incoming radiation.
%
%   INFO - string-description for the structure-IAM model, or just the ouput of func2str(IAM)
%
% See also: PVL_IAM_ASHRAE, PVL_IAM_PHYSICAL, GUIIRRTRANS, PVL_IAM_MARTINRUIZ, DIFFUSEIAM

    if nargin < 1 || isempty(IAM)
        if evalin('base','exist(''IAM'')')
        % Use IAM from base workspace only if checkIAM(IAM) runs without warnings
           IAM = evalin('base','IAM');
           warning_resetter = naptime({'checkIAM:defaults','checkIAM:weird',...
                        'completestruct:BnotinA','completestruct:BnotinA'},'error'); %#ok<NASGU>
           try
              [f,info,IAM] = checkIAM(IAM);
              return;
           catch
               warning('checkIAM:base',['IAM in base workspace failed to run smoothly, ',...
                        'falling back to SimOptions.IAM']);
           end
           clear lastwill
        end
        IAM = getSimOption('IAM'); 
    end
    if isequal(IAM,1) || ischar(IAM) && strcmpi(IAM,'none'), IAM = struct('model','none'); end
    
    % n = 1.526 for glass, L = 0.002 m, K = 4·1/m clear glass
    DEF.refraction = struct('model','refraction',...
        'refractionIndex',1.526,'glazingThickness',0.002,'extintionCoeff',4.0);
    
    % a_r = 0.17 for a clear mSi module, 0.20 with dust
    DEF.martinruiz = struct('model','martinruiz','ar',0.17,'maxloss',[]);
    
    DEF.hcpv = struct('model','hcpv','maxIncAngle',5.0);
    DEF.ashrae = struct('model','ashrae','b',0.05,'maxIncAngle',[]);
    DEF.none = struct('model','none');
    
    ALIAS = {'none',{'none','ideal','uniform','1'};
             'ashrae',{'ashrae'};
             'refraction',{'refraction','physical'};
             'martinruiz',{'martin-ruiz','martinruiz'};
             'hcpv',{'cpv','hcpv'}};

    [opt,varargin] = getflagoptions(varargin,{'-plot'});

    if ischar(IAM)
        IAM = parselist(IAM,ALIAS,'IAM model','');
        % IAM = DEF.(IAM);
        % [f,info] = checkIAM(IAM);
        % warning('checkIAM:defaults',['Using default IAM parameters for ' info]);
        
        IAM = getpairedoptions([{'model',IAM},varargin],fieldnames(DEF.(IAM)),'restchk');
    else
        assert(isempty(varargin),'Unrecognized arguments');
    end
        
    scalarposfld = @(s,fld) isfield(s,fld) && isscalar(s.(fld)) && s.(fld) > 0;
    if isstruct(IAM)
        assert(isscalar(IAM) && isfield(IAM,'model'),'Invalid structure');
        IAM.model = parselist(IAM.model,ALIAS,'IAM model','');
        DEF = DEF.(IAM.model);
        
        switch IAM.model
        case 'none'
            info = sprintf('No IAM correction!');
            f = @(x) 1*(x < 90);
            
        case 'ashrae'
            if scalarposfld(IAM,'maxIncAngle') 
                if ~scalarposfld(IAM,'b')
                    IAM.b = 1/(secd(IAM.maxIncAngle)-1);
                else
                    DEF.maxIncAngle = asecd(1 + 1./IAM.b);
                    if abs(DEF.maxIncAngle - IAM.maxIncAngle) > 0.1
                        warning('Replacing inconsistent MaxIncAngle %0.1f -> %0.1f',...
                                IAM.maxIncAngle,DEF.maxIncAngle);
                    end
                end
            end
            IAM = completestruct(IAM,DEF,'valid',@isscalar,'warning','A<B');
            IAM.maxIncAngle = asecd(1 + 1./IAM.b);

            b = IAM.b;
            f = @(x) pvlmod_iam_ashrae(b,min(90,abs(x)));
            info = sprintf('ASHRAE IAM model, b = %0.3f, th_max = %0.1f°',b,IAM.maxIncAngle);

        case 'refraction'
            IAM = completestruct(IAM,DEF,'valid',@isscalar,'warning','A<B');
            L = IAM.glazingThickness;
            K = IAM.extintionCoeff;
            n = IAM.refractionIndex;
            f = @(x) reshape(pvl_iam_physical(K,L,n,min(90,abs(x(:)))),size(x));
            info = sprintf('Physical-refraction IAM: K = %0.1f/m, L = %0.1fmm, n = %0.2f',K,L*1000,n);

        case 'martinruiz'
            if scalarposfld(IAM,'maxloss') 
                if ~scalarposfld(IAM,'ar')
                    IAM.ar = martinruizfit(IAM.maxloss);
                else
                    DEF.maxloss = martinruizmaxloss(IAM.ar);
                    if abs(DEF.maxloss - IAM.maxloss) > 1e-3
                        warning('Replacing inconsistent maxloss %0.1f%% -> %0.1f%%',...
                                IAM.maxloss*100,DEF.maxloss*100);
                    end
                end
            end
            IAM = completestruct(IAM,DEF,'valid',@isscalar,'warning','A<B');
            IAM.maxloss = martinruizmaxloss(IAM.ar);
            
            ar = IAM.ar;
            f = @(x) reshape(pvl_iam_martinruiz(ar,min(90,abs(x(:)))),size(x));
            info = sprintf('Martin-Ruiz IAM model: a_r = %0.2f (max-loss %0.1f%%)',ar,...
                IAM.maxloss*100);

        case 'hcpv'
            IAM = completestruct(IAM,DEF,'warning','A<B');
            thmax = IAM.maxIncAngle;
            f = @(x) 1*(x < thmax);
            info = sprintf('Threshold IAM model, th_max = %0.2f°',thmax);
        end
        
    elseif isa(IAM,'function_handle')
        
        f = IAM;
        try info = func2str(f); catch, info = 'Unknown function'; end
    else
        error('Expecting structure or function-handle');
    end
    
    msg = '';
    try 
        v = f(-180:180);
        validateattributes(v,{'numeric'},{'size',[1,361],'real','nonnegative','<',1.5});
    catch
    % see if clipping at 90° fixes the issue
        g = f;
        f = @(x) g(min(abs(x),90));
    end
    
    try 
        v = f(0:180);
        validateattributes(v,{'numeric'},{'size',[1,181],'real','nonnegative','<',1.5},'','f(0:90)');
    catch ERR
        msg = getReport(ERR);
    end
    assert(isempty(msg),'Bad function-handle: %s',msg);
    
    warnmsg = {};
    if f(0) ~= 1, warnmsg{1} = sprintf('f(0) = %0.2f %c 1',f(0),char('<'+2*(f(0)>1))); end
    if f(90) ~= 0, warnmsg{end+1} = 'f(90) ~= 0'; end
    if any(v < 0), warnmsg{end+1} = 'f(x) < 0 for some x in [0,180]'; end
    if any(diff(v) > 0), warnmsg{end+1} = 'f(x) is not monotonic descending'; end
    if ~isempty(warnmsg)
        warnmsg = sprintf('Unexpected IAM function: %s',strjoin(warnmsg,', '));
        if nargout < 4, warning('checkIAM:weird',warnmsg); end  %#ok<SPWRN>
    end

    if opt.plot, plotIAM(f,info); end
end

function plotIAM(fIAM,info)
% Plot beam- and integrated diffuse functions for fIAM

    [iso,alb,hb] = diffuseIAM(fIAM);
    
    GUIfigure('iam','Incident Angle Modifier'); clf(); hold on;
    tt = 0:90;
    plot(tt,fIAM(tt));
    plot(90-tt,iso(tt));
    plot(90-tt,alb(tt));
    plot(90-tt,hb(tt));
    axis([0,90,0.8,1]);
    legend({'Beam/CS','Isotropic Sky','Isotropic Albedo','Horizon Line'});
    legend box off
    xlabel('Beam Incidence Angle / 90° - Surface Tilt (degrees)');
    ylabel('IAM Factor');
    grid on
    title(info,'interpreter','none');
end

function [a,c] = martinruizfit(max_loss)
% [A,C] = MARTINRUIZFIT(L) - find the coefficient A for which the maximum energy loss due to IAM
%   is bound by L, that is:
%
%       max{(1-PVL_IAM_MARTINRUIZ(A,x))·cos(x)} = L  at x = C
%
% By construction, the maximum loss factor must be less than 0.25 (corresponding to A = Inf).
% For small values of A, the max. loss is approx. ~ A/exp(1) and takes place at ~80°

    validateattributes(max_loss,{'numeric'},{'scalar','positive','real','<',0.25},'','max_loss');

    a = fzero(@(x) martinruizmaxloss(x)-max_loss,exp(1)*max_loss);

    b = -1/a;
    x = (lambertwlog(b+1)-1)./b;
    c = acosd(x);
end

function y = martinruizmaxloss(a)

    if a < 0.05
        y = a./exp(1);
    else
        b = -1/a;
        x = (lambertwlog(b+1)-1)./b;
        y = -(b.*exp(b).*x.^2)./((x.*b+1).*(1-exp(b)));
    end
end