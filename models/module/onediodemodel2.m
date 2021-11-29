function [V,dVdI,d2VdI2] = onediodemodel2(P,Ix,Tol,MaxIter)
% Vx = ONEDIODEMODEL2(PARAMS,Ix,TOL,MAXITER) - Calculates the solution Vx(Ix), to the 
%   implicit one-diode-model equation, including an optional recombination- and breakdown current
%   terms. See details about notation and parameter input in ONEDIODEMODEL.
% PARAMS: structure of parameters, usually generated by TRANSLATEODM.
% Ix: numeric array [A], any size (results will have the same size).
% TOL: maximum absolute tolerance of the result. Default is Iph·SimOptions.RelTol, or zero*.
%      NOTE: the minimum achievable tolerance is bound by NEPS·eps(Vx) (see below).
% MAXITER: optionally sets iteration limit for final Newton-Raphson iteration.
% Vx: size(Ix) solution array.
%
% See also: ONEDIODEMODEL, TRANSLATEODM, ONEDIODEMPP, ONEDIODEPWLAPPROX, CHECKODM

    persistent OPT
    if isempty(OPT), OPT.MaxIter = getSimOption('MaxIter'); end
    
    narginchk(2,4);    
    assert(isnumeric(Ix) && isreal(Ix),'onediodemodel2:Ix','Ix must be numeric and real');
    
    if nargin < 4, MaxIter = OPT.MaxIter; end
    if nargin < 3 || isempty(Tol), Tol = []; end
    [~,tolv,toli] = parsetolerance(Tol);
    
    i0 = P.Iph*(1-P.di2mutau/P.Vbi);  % I @ Vj = 0
	fwdb = (Ix <= i0);                % Vj > 0
    revb = ~fwdb & isfinite(Ix);
    
    % Solve ODM in forward bias, getting dIdV, d2IdV2 only if asked
    out = cell(1,max(nargout(),1));
    [out{:}] = forwardODM2(P,Ix(fwdb),tolv,toli,MaxIter);
    V = NaN(size(Ix));
    V(fwdb) = out{1};
    if nargout > 1
        dVdI = NaN(size(Ix));
        d2VdI2 = NaN(size(Ix));
        dVdI(fwdb) = 1./out{2};                 % dx/dy = 1/(dx/dy)
        d2VdI2(fwdb) = -out{3}.*dVdI(fwdb).^3;  % d²x/dy² = -d²y/dx²·(dx/dy)³
    end
    
    % Do the same for reverse bias
    if ~all(fwdb,'all')
        if all(isfield(P,{'avalanche_a','avalanche_m','Vbr'}))
            [out{:}] = reverse_bishop2(P,Ix(revb),tolv,toli);
        elseif isfield(P,'Brev')
            [out{:}] = reverse_pvsyst2(P,Ix(revb),tolv,MaxIter);
        elseif all(isfield(P,{'Vbr','Vbi'}))
            [out{:}] = reverse_alonso2(P,Ix(revb),tolv,toli);
        else
            [out{:}] = forwardODM2(P,Ix(revb),tolv,MaxIter);
            warning('odm:revmodel','No reverse-bias parameters!');
        end
        V(revb) = out{1};
        if nargout > 1
            dVdI(revb) = 1./out{2};
            d2VdI2(revb) = -out{3}.*dVdI(revb).^3;
        end
    end
end

function [v,di,d2i] = forwardODM2(P,i,tolv,toli,MaxIter)
% Solve the 7-parameter diode model

%     F = ones(size(Ix));
    i0 = P.Iph*(1-P.di2mutau/P.Vbi);

    if (P.di2mutau > 0) && P.Iph > 0
    % Breakdown voltage Vbi effectively limits the solution space to [-Rs·i0,Vbi)
        err = @(v) P.Iph-exp(log(P.Io)+(v+i*P.Rs)/P.nVth)+P.Io -(v+i*P.Rs)/P.Rsh - ...
                   P.Iph*P.di2mutau./max(eps(0),P.Vbi-v-P.Rs*i) - i;
        
        v = bisection(err,-i0*P.Rs,P.Vbi-P.Rs*i,tolv,toli);
%         F(~revb) = 0.2;
    else
    % Solve the 5-parameter ODM for positive-biased points
        v = fiveparameterODM2(i,P.Iph,P.Io,P.Rs,P.nVth,P.Rsh,tolv,MaxIter);
    end
    
% 	% Refine recombination/reverse-voltage currents (solve iteratively)
%     if any(revb) || (P.di2mutau > 0)
%         for iter = 1:MaxIter+1
% 			oldVx = Vx;
%             Irec = P.Iph*P.di2mutau./(P.Vbi-Vx-P.Rs*Ix);
%             Irev = revb.*P.Brev.*(Vx+P.Rs*Ix).^2;
% 			newIx = fiveparameterODM(Vx,P.Iph-Irec+Irev,P.Io,P.Rs,P.nVth,P.Rsh,toli,MaxIter);
%             err = newIx - Ix;
%             dIx = dIdV(P,Vx,newIx,revb);
%             Vx = Vx - F.*err./dIx;         
%             if all(abs(Vx-oldVx) <= tolv(Vx) | abs(err) <= toli(Ix)), break; end %(**)
%         end
%         if iter > MaxIter
%             error('ODM2 could not converge to a solution'); 
%         end
%     end
% 
%     function dIx = dIdV(P,Vx,Ix,revb)
%         % Explicit 1st and 2nd derivatives
%         kc = P.Io*P.Rsh/P.nVth*exp((Vx + P.Rs*Ix)/P.nVth);
%         kd = (P.Iph*P.Rsh*P.di2mutau)./(Vx - P.Vbi + P.Rs*Ix).^2;
%         ke = kc + kd + 1;
%         dIx = -ke./(P.Rsh + P.Rs*ke);
%         % d2Ix = (1+P.Rs*dIx).^2.*(2*kd./(Vx - P.Vbi + P.Rs*Ix)-kc/P.nVth)./(P.Rsh + P.Rs*ke);
% 
%         kr = 2*P.Brev*(Vx+P.Rs*Ix).*revb;
%         dIx = (dIx + kr)./(1 - P.Rs*kr);
%         %d2Ix(revb) = (d2Ix(revb) + 2*P.Brev*(1-P.Rs*dIx(revb)).^2)./(1 - P.Rs*kr);
%     end

    % Explicit 1st and 2nd derivatives
    if nargout > 1
        [di,d2i] = ODMderivatives(P,v + P.Rs*i);
    end
end

function V = fiveparameterODM2(I,Iph,Io,Rs,nVth,Rsh,tolv,MaxIter)
% Finds the solution to the basic ODM equation: i = Iph-io*(exp((vx+i*Rs)/nVth)-1)-(vx+i*Rs)/Rsh
% using the Lambert-W function... or iteratively, when this fails

    % Solve using Lambert-W function (real,positive branch)
    %   th =ka*exp(ka*Idf/Io), where ka = Io*Rsh/nVth, Idf = Iph-Ix
    %   Vx = Rsh*Idf-Rs*Ix-nVth*W(th)

    V = zeros(size(I));

    ka = Io*Rsh/nVth;
    Idf = Iph-I;

    % Module is Reverse-Biased from the point where I > Iph
    rb = (Idf <=0);
    V(rb) = Iph-Idf(rb)*Rs;

    logth = log(ka)+(ka*Idf(~rb)/Io); % th can be huge, so get ln(th) first
    W = lambertwlog(logth);
    V(~rb) = Rsh*Idf(~rb)-Rs*I(~rb)-nVth*W;

    % lambertwlog ensures precision within O(eps(W)) for W, but this doesn't imply sufficient
    % precission in Vx. Perform a few Newton-Raphson iterations directly over Vx, if necessary:
    % NOTE: use of exp(log(Io)+ ...) instead of Io·exp(...) helps reduce rounding error
    
    f = @(Vx) Rsh*(Iph-exp(log(Io)+(Vx+I*Rs)/nVth)+Io)-I*(Rs+Rsh)-Vx;
    df = @(Vx) -Rsh/nVth*exp(log(Io)+(Vx+I*Rs)/nVth) - 1;
    if any(abs(f(V)) > tolv(V)) %(**)
        V = quickfzero(f,df,V,tolv,MaxIter);
    end
end

function [v,di,d2i] = reverse_bishop2(P,i,tolv,toli)
% Bishop (1988) exponential model, solved by bisection, knowing that junction voltage must lie
% between v + i(vj = 0)·Rs and 0

    err = @(vj) P.Iph*(1-P.di2mutau./(P.Vbi-vj)) - vj/P.Rsh.*(1 + P.avalanche_a*(1-vj/P.Vbr).^(-P.avalanche_m))-i;
    vj = bisection(err,P.Vbr+i*P.Rs,0,tolv,toli);
    v = vj - i*P.Rs;
    
    if nargout > 1
        K = 1./(P.Rsh*(P.Vbr - vj)).*P.avalanche_a.*(1-vj/P.Vbr).^(-P.avalanche_m);
        djA = K.*(P.Vbr + (P.avalanche_m-1)*vj);
        d2jA = K.*P.avalanche_m.*(2*P.Vbr + (P.avalanche_m-1)*vj)./(P.Vbr - vj);
        [di,d2i] = ODMderivatives(P,vj,djA,d2jA);
    end
end

function [v,di,d2i] = reverse_alonso2(P,i,tolv,toli)
% Alonso-García & Ruiz (2006) model, explicit in terms of V
    
    Be = 3.0;
    i0 = P.Iph*(1-P.di2mutau/P.Vbi);    % I @ Vj = 0
    m0 = ODMderivatives(P,0);           % dI/dV @ Vj = 0
    
    % iN(v) = Linear approx. to I (w/o Ia) for V < -IRs
    % NOTE: i0 is adjusted by (1-Ke0) so that @ Vj = 0, iN/(1-Ke) = i0 
    Ke0 = exp(Be*(1-sqrt((P.Vbi-P.Vbr)./(P.Vbi+i0*P.Rs))));
    iN = @(v) (1-Ke0)*i0 + m0.*(v+i0*P.Rs);  
    
    % Search over v, min. voltage bound by the condition iN < i
    err = @(v) P.Vbr + (P.Vbi - P.Vbr)*(1-(1-log(1-iN(v)./i)/Be).^(-2))-v;
    vm = max(P.Vbr,(i-i0*(1-Ke0))/m0-i0*P.Rs); 
    v = bisection(err,0,vm,tolv,toli);
        
    if nargout > 1
        F = sqrt((P.Vbi-P.Vbr)./(P.Vbi-v));
        Ke = exp(Be*(1-F));
        k = Be*Ke.*F./((P.Vbi-v).*(1-Ke).^2);
        iN = iN(v);
        di = m0./(1-Ke) - 0.5*iN.*k;
        d2i = k.*(iN./(4*(P.Vbi-v)).*(Be*F.*(1+Ke)./(1-Ke)-3)-m0);
    end
end

function [v,di,d2i] = reverse_pvsyst2(P,i,tolv,MaxIter)
% PVsyst quadratic model: solution to I = Iph + ... + Brev(V+RsI)²


    i0 = P.Iph*(1-P.di2mutau/P.Vbi);    % I @ Vj = 0
    m0 = ODMderivatives(P,0);           % dI/dV @ Vj = 0

    % Solution to I = i0 + m0(V+RsI) + Brev(V+RsI)²
    v = -P.Rs*i - m0/(2*P.Brev)*(1-sqrt(1+4*P.Brev/m0*(i-i0).*(1/m0+P.Rs)));

    if P.di2mutau > 0
    % Fix deviations due to recombination current, if necessary
        f = @(v) P.Iph*(1-P.di2mutau./(P.Vbi-v-i*P.Rs))-(v+i*P.Rs)/P.Rsh + P.Brev*(v+i*P.Rs).^2-i;
        df = @(v) -P.Iph*P.di2mutau./(P.Vbi-v-i*P.Rs).^2 -1/P.Rsh + 2*P.Brev*(v+i*P.Rs);
        v = quickfzero(f,df,v,tolv,MaxIter);
    end
  
    if nargout > 1
        vj = v+i*P.Rs;
        djA = -2*P.Brev*vj;
        d2jA = -2*P.Brev;
        [di,d2i] = ODMderivatives(P,vj,djA,d2jA);
    end
end