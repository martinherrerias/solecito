function varargout = onediodepwlapprox(P,Tol,Lim,varargin)
% [V,I] = ONEDIODEPWLAPPROX(PARS,TOL,LIM) - Selects a set of points V,I over the ODM curve 
%   defined by parameters PARS by recursively dividing interval LIM until a piece-wise-linear
%   approximation within all intervals is reasonable within TOL. The starting interval 
%   [a,b] is provided by LIM (see below).
%
%   NOTES: The algorithm works on an explicit representation of the ODM equations using junction
%   voltage V + IRs as an auxiliary variable (tolerance in X is also not strictly tolerance in V).
%   If the inflection point V = -IRs is included in LIM, the algorithm makes sure to split
%   the problem into two convex-approximation sub-problems.
%
%   [V,I] = ONEDIODEPWLAPPROX(PARS,TOL,LIM,VARARGIN) - Pass additional arguments to PWLAPPROX,
%   such as partition criteria, error metric, and/or bias-correction factor.
%
% OBJ = ONEDIODEPWLAPPROX(PARS,TOL,LIM,..) - return an MDPWL object instead of V,I vectors.
%
%   PARS: ODM parameter set - vector [P.Iph P.Io P.Rs P.nVth P.Rsh] and optionally [... P.di2mutau P.Vbi P.Brev]
%         or structure with those field names, usually generated by TRANSLATEODM
%   TOL: (incomplete) vector of approximation tolerances. See PARSETOLERANCE.
%   LIM: 2-vector defining voltage limits for the interval. Default is [0,Voc], or...
%        4-vector of [Vmin,Vmax,Imin,Imax], equivalent to [max(Vmin,V(Imax)),min(Vmax,V(Imin))]
%        i.e. restrictions are additive, LIM defines a 'crop window'.
%
% Example: [v,i] = ONEDIODEPWLAPPROX([10 1e-8 0.3 1.8 400],1e-3); figure(); plot(v,i,'ro');
%          [v,i] = ONEDIODEPWLAPPROX([10 1e-8 0.3 1.8 400],1e-3); figure(); plot(v,i,'ro');
%
% See also: PWLAPPROX, ONEDIODEMODEL, TRANSLATEODM, MODULEINTERPOLANT, MDPWL

    narginchk(1,7);
    assert(isstruct(P) && isscalar(P),'Expecting scalar structure');
    % parsestruct(P,{'Iph','Io','Rs','nVth','Rsh'},'opt',{'di2mutau','Vbi','Brev'},'-n','-r','-s');
    if nargin < 2 || isempty(Tol), Tol = []; end
    
    % Set interval as [0, Voc,-Inf,Inf], when omitted
    if nargin < 3, Lim = []; end
    
    % varargin = {method,metric,slope,bias_K}
    varargin(end+1:4) = {[]};
    if isempty(varargin{3}), varargin{3} = '2nd'; end
            
    if isempty(Tol) || isempty(Lim)
        [~,Vmp,Imp,Voc,Isc] = onediodeMPP(P,0);
    end
    
    if isempty(Tol)
        Tol = parsetolerance(Tol);  % just to get Tol(3) = RelTol
        Tol = parsetolerance([Vmp/10,Imp/10,1]*Tol(3));
    else
        Tol = parsetolerance(Tol);
    end

    switch numel(Lim)
        case 0
            LimV = [0 Voc];
            LimI = [Isc 0];
        case 2
            LimV = sort(Lim(:))';
            LimI = onediodemodel(P,LimV,0); % ($)
        case 4
            Lim(~isfinite(Lim)) = NaN;
            LimV = [Lim(1:2),onediodemodel2(P,Lim(3:4),0)];
            LimV = [max(LimV([1,4]),[],'omitnan'),min(LimV([2,3]),[],'omitnan')];
            LimI = onediodemodel(P,LimV,0);
        otherwise, LimV = [0,0]; % make it crash below
    end
    assert(LimV(2) > LimV(1),'onediodepwlapprox:lim','Bad/empty limits');
	
    % Move to auxiliary variable z = V + IRs (junction voltage)
    LimVj = (LimV + P.Rs*LimI);
    
    % Pick reverse-bias model
    if all(isfield(P,{'avalanche_a','avalanche_m','Vbr'})), mdl = 'b'; % Bishop
        elseif isfield(P,'Brev'), mdl = 'p'; % PVsyst
        elseif all(isfield(P,{'Vbr','Vbi'})), mdl = 'a'; % Alonso
        else, mdl = 'n'; % none
    end
        
    % Get PWL approximation of I(z) in interval LimZ...
    
    if LimVj(1) >= 0 || mdl == 'n'
    % Diode is forward biased in complete interval (or there's no reverse-bias model)
        [vj,I] = pwlapprox(@explicitODM_fwd,LimVj(1),LimVj(2),Tol,varargin{:});
        V = vj - I*P.Rs;
    else
        if LimVj(2) <= 0
        % Diode is reverse biased in complete interval
            switch mdl
            case 'b'
                [vj,I] = pwlapprox(@reverse_bishop_vj,LimVj(1),LimVj(2),Tol,varargin{:});
            case 'p'
                % NOTE: PVsyst reverse model generates an inflection point at ...
                vj0 = P.Vbi - (P.Iph*P.di2mutau/P.Brev)^(1/3);
                if vj0 < LimVj(2) && vj0 > LimVj(1)
                    [Vj{1},i{1}] = pwlapprox(@reverse_pvsyst_vj,LimVj(1),vj0,Tol,varargin{:});
                    [Vj{2},i{2}] = pwlapprox(@reverse_pvsyst_vj,vj0,LimVj(2),Tol,varargin{:});
                    vj = [Vj{1}(1:end-1),Vj{2}];
                    I = [i{1}(1:end-1),i{2}];
                else
                    [vj,I] = pwlapprox(@reverse_pvsyst_vj,LimVj(1),LimVj(2),Tol,varargin{:});
                end
            case 'a'
                [V,I] = pwlapprox(@reverse_alonso,LimV(1),LimV(2),Tol,varargin{:});
            end
            if mdl ~= 'a', V = vj - I*P.Rs; end
        else
        % Interval contains vj = 0            
            switch mdl
                case 'b'
                    [vj,i{1}] = pwlapprox(@reverse_bishop_vj,LimVj(1),0,Tol,varargin{:});
                case 'p'
                    % NOTE: PVsyst reverse model generates an inflection point at ...
                    vj0 = P.Vbi - (P.Iph*P.di2mutau/P.Brev)^(1/3);
                    if vj0 < 0 && vj0 > LimVj(1)
                        [Vj{1},i{1}] = pwlapprox(@reverse_pvsyst_vj,LimVj(1),vj0,Tol,varargin{:});
                        [Vj{2},i{2}] = pwlapprox(@reverse_pvsyst_vj,vj0,0,Tol,varargin{:});
                        vj = [Vj{1}(1:end-1),Vj{2}];
                        i{1} = [i{1}(1:end-1),i{2}];
                    else
                        [vj,i{1}] = pwlapprox(@reverse_pvsyst_vj,LimVj(1),0,Tol,varargin{:});
                    end
                case 'a'
                    v0 = -P.Rs*P.Iph*(1-P.di2mutau/P.Vbi);    % V @ Vj = 0
                    [v{1},i{1}] = pwlapprox(@reverse_alonso,LimV(1),v0,Tol,varargin{:});
            end
            if mdl ~= 'a', v{1} = vj - i{1}*P.Rs; end
            
            [vj,i{2}] = pwlapprox(@explicitODM_fwd,0,LimVj(2),Tol,varargin{:});
            v{2} = vj - i{2}*P.Rs;
        
            V = [v{1}(1:end-1),v{2}];
            I = [i{1}(1:end-1),i{2}];
        end
    end

    switch nargout
        case {2}, varargout = {V,I};
        case {0,1}, varargout = {mdpwl(V,I,eps(0))};
        otherwise, error('Too many output arguments')
    end
    
    function [I,di,d2i] = explicitODM_fwd(vj)
    % Explicit solution to the One-Diode-Model equation (forward bias), with recombination current, 
    % using as input the junction voltage Vj = V + I·P.Rs >= 0

        I = P.Iph*(1-P.di2mutau./(P.Vbi-vj))-P.Io*(exp(vj/P.nVth)-1)-vj/P.Rsh;

        % Explicit 1st and 2nd derivatives
        if nargout > 1
            [di,d2i] = ODMderivatives_vj(P,vj);
        end
    end

    function [i,di,d2i] = reverse_bishop_vj(vj)
    % Bishop (1988) exponential model, solved by bisection, knowing that junction voltage must lie
    % between v + i(vj = 0)·Rs and 0

        K = 1./P.Rsh.*P.avalanche_a.*(1-vj/P.Vbr).^(-P.avalanche_m);
        i = P.Iph*(1-P.di2mutau./(P.Vbi-vj)) - vj/P.Rsh - vj.*K;

        if nargout > 1
            K = K/(P.Vbr - vj);
            djA = K.*(P.Vbr + (P.avalanche_m-1)*vj);
            d2jA = K.*P.avalanche_m.*(2*P.Vbr + (P.avalanche_m-1)*vj)./(P.Vbr - vj);
            [di,d2i] = ODMderivatives_vj(P,vj,djA,d2jA);
        end
    end

    function [i,di,d2i] = reverse_pvsyst_vj(vj)
    % PVsyst quadratic model: I = Iph + ... + Brev(V+RsI)²

        i = P.Iph*(1-P.di2mutau./(P.Vbi-vj)) - vj/P.Rsh + P.Brev*vj.^2;

        if nargout > 1
            djA = -2*P.Brev*vj ;
            d2jA = -2*P.Brev;
            [di,d2i] = ODMderivatives_vj(P,vj,djA,d2jA);
        end
    end
 
    function [i,di,d2i] = reverse_alonso(v)
    % Alonso-García & Ruiz (2006) model, explicit in terms of V

        v(v < P.Vbr) = NaN;
        Be = 3.0;

        i0 = P.Iph*(1-P.di2mutau/P.Vbi);    % I @ Vj = 0
        m0 = ODMderivatives_vj(P,0);        % dI/dV @ Vj = 0

        % iN = Linear approx. to I (w/o Ia) for V < -IRs
        % i0 is adjusted by (1-Ke0) so that @ Vj = 0, iN/(1-Ke) = i0 
        Ke0 = exp(Be*(1-sqrt((P.Vbi-P.Vbr)./(P.Vbi+i0*P.Rs))));
        iN = (1-Ke0)*i0 + m0.*(v+i0*P.Rs);  

        F = sqrt((P.Vbi-P.Vbr)./(P.Vbi-v));
        Ke = exp(Be*(1-F));
        i = iN./(1-Ke);

        if nargout > 1
            k = Be*Ke.*F./((P.Vbi-v).*(1-Ke).^2);
            di = m0./(1-Ke) - 0.5*iN.*k;
            d2i = k.*(iN./(4*(P.Vbi-v)).*(Be*F.*(1+Ke)./(1-Ke)-3)-m0);
        end
    end
end

function [djI,d2jI] = ODMderivatives_vj(P,Vj,djA,d2jA)
% [dIdV,d2IdV2] = ODMDERIVATIVES(P,Vj,djA,d2jA) - Calculates first and second derivatives dI/dVj
%   and d²I/dVj² for the ODM curve at junction voltages Vj, optionally considering an 
%   additional current term Ia (e.g. avalanche breakdown on reverse bias), or more specifically, 
%   its derivatives with respect to Vj: djA = dIa/dVj, and d2jA = d²Ia/dVj²

    if nargin < 3, djA = 0; end
    if nargin < 4, d2jA = 0; end
        
    djD = P.Io/P.nVth*exp(Vj/P.nVth);                % d(Idiode)/dVj
    djR = (P.Iph*P.di2mutau)./(P.Vbi-Vj).^2;         % d(Irec)/dVj
    djI = -(djD + djR + 1/P.Rsh + djA);              % dI/dVj    
    d2jI = -(2*djR./(P.Vbi-Vj)+djD/P.nVth + d2jA);   % d²I/dVj²
end


