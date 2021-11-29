function [Pmp,Vmp,Imp,Voc,Isc,c] = onediodeMPP(Params,tol,Vmp_ref,MaxIter)
% [Pmp,Vmp,Imp,Voc,Isc,c] = ONEDIODEMPP(PARAMS,TOL,MAXITER)  - Find the MPP of a one-diode model 
%   function defined by PARAMS, within precision TOL. Return also Voc, Isc, and c* as byproducts.
%
%  PARAMS: (array of) structure(s) of 5-8 ODM parameters (see ONEDIODEMODEL for details)
%  TOL: (incomplete) unified tolerance vector. See PARSETOLERANCE
%  MAXITER: optional limit of iterations for Newton-Raphson search.
%
%  [Pmp,Vmp,Imp,Voc,Isc,c] - size PARAMS arrays of MPP power, voltage, and current; open-circuit
%  voltage, short-circuit current, and curvature parameter c*, respectively.
%
%  (*) Parameter c = Vmp·(d²I/dV²)/(dI/dV) = -Vmp²/Imp·(d²I/dV²) @ MPP, can be used to get a 2nd 
%  order approximation of the IV curve in the vecinity of Vmp, namely:
%
%       I/Imp = a + b·exp(c·V/Vmp), with a = (1+1/c), and b = exp(-c)/c
%
%  c is useful to estimate the impact of small mismatch losses, proportional to (c+2)/c. See
%  L. L. Bucciarelli, “Power loss in photovoltaic arrays due to mismatch in cell characteristics” 
%  Solar Energy, vol. 23, no. 4, pp. 277–288, Jan. 1979, doi: 10.1016/0038-092X(79)90121-X.
%
% See also: ONEDIODEMODEL, TRANSLATEODM, PARSETOLERANCE
    
    persistent OPT
    if isempty(OPT)
        OPT.MaxIter = getSimOption('MaxIter');
        OPT.NTOL = 16;  % tolerance refinement factor for calls to ONEDIODEMODEL
    end
    
    if nargin < 2, tol = 0; end
    tol = parsetolerance(tol/OPT.NTOL);             % make tol = tol/OPT.NTOL for ONEDIODEMODEL calls
    [~,tolv,toli] = parsetolerance(tol*OPT.NTOL);   % ensure MPP tolerance is actually OPT.NTOL*tol
    
    if nargin < 4, MaxIter = OPT.MaxIter; end
    if nargin < 3, Vmp_ref = []; end
    if ~isempty(Vmp_ref)
        assert(numel(Vmp_ref) == numel(Params),'If provided, Vmp_ref must have the size of PARAMS'); 
    end
    
    if numel(Params) > 1
        if isempty(Vmp_ref)
            [Pmp,Vmp,Imp,Voc,Isc,c] = arrayfun(@(p) onediodeMPP(p,tol,Vmp_ref,MaxIter),Params);
        else
            [Pmp,Vmp,Imp,Voc,Isc,c] = arrayfun(@(p,v) onediodeMPP(p,tol,v,MaxIter),Params,Vmp_ref);
        end
        return;
    end
    
    % Get Voc,Isc with increased precission
    Voc = onediodemodel2(Params,0,tol);
    Isc = onediodemodel(Params,0,tol);

    % There is an actual maximum at Vmp < 0 for Voc < 0 and Isc < 0, but who cares
    if Isc <= toli(0) || Voc <= tolv(0)
        Vmp = 0; Imp = Isc; Pmp = 0; Voc = 0; Isc = 0; c = 0;
        return
    end
    
    % Debug
    % figure
    % plot(0:Voc/99:Voc,onediodemodel(Params,0:Voc/99:Voc));
    % hold on
    
    if ~isempty(Vmp_ref) && Vmp_ref > 0 && Vmp_ref < Voc
        Vmp = Vmp_ref;
    else
    % NOTE: 0.81·Voc is the median for most modules, but convergence is faster fromt the right
        Vmp = Voc*0.84; % Seed value 
    end
    
    [Imp,di,d2i] = onediodemodel(Params,Vmp,tol);
    for iter = 0:MaxIter 
        oldv = Vmp; %oldp = Pmp;
          
        if iter < 2
        % MPP for the approximation i/i0 = a - b·exp(cv/v_ref) [inspired by Bucciarelli, 1979]
        % NOTE: once LAMBERTWLOG is optimized/translated, it might be convenient to use this
        % approximation all the way. 
            c = Vmp*d2i/di;
            k = (1-d2i*Imp/di^2);
            Vmp = (lambertwlog(log(k)+c+1)-1)*Vmp/c;
        else
        % 2nd order Taylor expansion solution
            Vmp = Vmp-(Vmp*di+Imp)/(Vmp*d2i+2*di);
        end
        [Imp,di,d2i] = onediodemodel(Params,Vmp,tol);
        
        if abs(oldv-Vmp) <= tolv(Vmp) && abs((oldv-Vmp)*di) <= toli(Imp), break; end
    end
    if iter >= MaxIter
        % Vmp = NaN; Imp = NaN; Pmp = NaN;
        error('onediodeMPP:itermax','Maximum number of iterations reached');
    end
    Pmp = Vmp*Imp;
    c = Vmp*d2i/di;
end
