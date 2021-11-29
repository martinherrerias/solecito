function [Pmp,Vmp,Imp,Voc,Isc] = onediodeMPP(Params,tol,Vmp_ref,MaxIter)
% [Pmp,Vmp,Imp,Voc,Isc] = ONEDIODEMPP(PARAMS,TOL,MAXITER)
% Finds the MPP of a one-diode model function defined by Params, within precision TOL.
%
%  PARAMS: (array of) structure(s) of 5-8 ODM parameters (see ONEDIODEMODEL for details)
%  TOL: (incomplete) unified tolerance vector. See PARSETOLERANCE
%  MAXITER: optional limit of iterations for Newton-Raphson search.
%
%  [Pmp,Vmp,Imp,Voc,Isc] - size PARAMS arrays of MPP power, voltage, and current; open-circuit
%  voltage, and short-circuit current, respectively.
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
        Vmp = zeros(size(Params)); Imp = Vmp; Pmp = Vmp; Voc = Pmp; Isc = Pmp;
        for j = 1:numel(Params)
            if isempty(Vmp_ref), Vmp_refj = []; else, Vmp_refj = Vmp_ref(j); end
            [Pmp(j),Vmp(j),Imp(j),Voc(j),Isc(j)] = onediodeMPP(Params(j),tol,Vmp_refj,MaxIter);
        end
        return;
    end
    
    % Get Voc,Isc with increased precission
    Voc = onediodemodel2(Params,0,tol);
    Isc = onediodemodel(Params,0,tol);

    % There is an actual maximum at Vmp < 0 for Voc < 0 and Isc < 0, but who cares
    if Isc <= toli(0) || Voc <= tolv(0)
        Vmp = 0; Imp = Isc; Pmp = 0;
        return
    end
    
    % Debug
    % figure
    % plot(0:Voc/99:Voc,onediodemodel(Params,0:Voc/99:Voc));
    % hold on

    if ~isempty(Vmp_ref) && Vmp_ref > 0 && Vmp_ref < Voc
        Vmp = Vmp_ref;
    else
        Vmp = Voc*0.85;      % Seed value = 0.85·Voc
    end
    
    [Imp,di,d2i] = onediodemodel(Params,Vmp,tol);
    for iter = 0:MaxIter 
        oldv = Vmp; %oldp = Pmp;
        Vmp = Vmp-(Vmp*di+Imp)/(Vmp*d2i+2*di);
        [Imp,di,d2i] = onediodemodel(Params,Vmp,tol);
        if abs(oldv-Vmp) <= tolv(Vmp) && abs((oldv-Vmp)*di) <= toli(Imp), break; end
    end
    if iter >= MaxIter
        % Vmp = NaN; Imp = NaN; Pmp = NaN;
        error('onediodeMPP:itermax','Maximum number of iterations reached');
    end
    Pmp = Vmp*Imp;
end
