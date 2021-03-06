function [I,dIdV,d2IdV2] = onediodemodel(P,V,Tol,MaxIter)
% [I,dIdV,d2IdV2] = ONEDIODEMODEL(PARAMS,V,TOL,MAXITER) - Calculates the solution i(v) to the 
%   implicit one-diode-model equation, including an optional recombination-current term [1],[4]:
%
%       i = Iph-Io·(exp(vj/nVth)-1)-vj/Rsh - Iph·di2mutau/(Vbi-vj) - Ia(v,i)
%       vj = v + i·Rs
%
%   Where:
%     Iph = photocurrent [A]
%     Io = dark or diode reverse saturation current [A]
%     Rs = Series resistance [Ohm]
%     Rsh = Shunt resistance [Ohm]
%     nVth = equivalent thermal voltage [V], product of the diode ideality factor (n), 
%       number of cells in series (Ns), and cell thermal voltage.
%     di2mutau = recombination factor for amorphous cells [~1.4 V], ratio of di² (squared 
%       thickness of intrinsic layer), and effective drift-length (mu·tau)_eff. [1][4].
%     Vbi = Junction built-in voltage (all cells) [V], ~(0.9 V)·Ns, see [1],[4]
%     Ia(v,i) = breakdown current, at reverse-biased conditions (vj < 0). Calculated with one
%       of the following models (basd on parameter availability, with this precedence):
%
%       1. Bishop (1988) exponential model: Ia = Vj/Rsh·a·(1-Vj/Vbr)^(-m),
%           required parameters are avalanche-fraction a, exponent m, and breakdown voltage 
%           Vbr [V] named {avalanche_a, avalanche_m, Vbr} 
%       
%       2. PVsyst quadratic approximation: Ia = Brev·Vj²
%           required empirical reverse-bias coefficient Brev [A/V²]
%
%       3. Alonso-García & Ruiz (2006): i = i0(Vj)/(1-exp(Be·(1-sqrt((Vbi-Vbr)/(Vbi-V))))),
%           required parameters are built-in- and breakdown-voltages [V] {Vbi,Vbr}
%           NOTE: For this model, I(Vbr) = Inf, and I(V < Vbr) = NaN
%
%       4. No model. Ia = 0, with warning.
%
% INPUT:
%   PARAMS: structure with fields {Iph Io Rs nVth Rsh} and optionally {di2mutau, Vbi} to include 
%       a recombination term (otherwise di2mutau = 0). For reverse voltage characteristic, and
%       depending on the model to be used (see above) additional fields can include 
%       {avalanche_a, avalanche_m, Vbr, Vbi, Brev}. PARAMS is usually generated by TRANSLATEODM.
%   V: numeric vector of voltages, any size (results will have the same size).
%   TOL: (incomplete) unified tolerance vector. See PARSETOLERANCE
%       TOL - can be a vector of absolute tolerances [tolV,tolI]; a scalar relative tolerance, —   
%       in which case tolV = max(|V|·tol,eps(V)), and similarly for tolI; or it can be a  
%       3-vector [mintolV,mintolI,reltol]. In this case, tolV = max(mintolx,|x|·reltol,eps(x)),  
%       and similarly for y.
%   MAXITER: optionally sets iteration limit for final Newton-Raphson iteration.
%
% OUTPUT:
%   I, dIdV, d2IdV2: size(V) solution arrays: i, di/dx, and d²i/dx², respectively. The first and 
%                  second explicit derivatives are provided (e.g. for use in ONEDIODEMPP or
%                  ONEDIODEPWLAPPROX)
%
% [1] Mermoud, A., Lejeune, T., 2010. Performance assessment of a simulation model for PV modules 
%     of any available technology, in: 25th European PV Solar Energy Conference, Valencia.
% [2] Roger, J.A., Maguin, C., 1982. Photovoltaic solar panels simulation including dynamical 
%     thermal effects. Solar Energy 29, 245–256. https://doi.org/10.1016/0038-092X(82)90210-9
% [3] Merten, J., Asensi, J.M., Voz, C., Shah, A.V., Platz, R., Andreu, J., 1998. Improved 
%     equivalent circuit and analytical model for amorphous silicon solar cells and modules. 
%     Electron Devices, IEEE Transactions on 45, 423–429.
% [2] Bishop, J.W., 1988. Computer simulation of the effects of electrical mismatches in 
%   photovoltaic cell interconnection circuits. Solar cells 25, 73–89.
% [3] Alonso-García, M.C., Ruíz, J.M., 2006. Analysis and modelling the reverse characteristic of 
%   photovoltaic cells. Solar Energy Materials and Solar Cells 90, 1105–1120.
%
% See also: ONEDIODEMODEL2, TRANSLATEODM

    persistent OPT
    if isempty(OPT), OPT.MaxIter = getSimOption('MaxIter'); end

    narginchk(2,4);
    assert(isnumeric(V) && isreal(V),'onediodemodel2:I','I must be numeric and real');
        
    if nargin < 4, MaxIter = OPT.MaxIter; end
    if nargin < 3 || isempty(Tol), Tol = 0; end
    [~,~,toli] = parsetolerance(Tol);
    	
    i0 = P.Iph*(1-P.di2mutau/P.Vbi);  % I @ Vj = 0
	fwdb = (V >= -P.Rs*i0);           % Vj > 0
    revb = ~fwdb & isfinite(V);
    
    % Solve ODM in forward bias, getting dIdV, d2IdV2 only if asked
    out = cell(1,max(nargout(),1));
    [out{:}] = forwardODM(P,V(fwdb),toli,MaxIter);
    I = NaN(size(V));
    I(fwdb) = out{1};
    if nargout > 1
        dIdV = NaN(size(V));
        d2IdV2 = NaN(size(V));
        dIdV(fwdb) = out{2};
        d2IdV2(fwdb) = out{3};
    end
    
    % Do the same for reverse bias
    if ~all(fwdb,'all')
        if all(isfield(P,{'avalanche_a','avalanche_m','Vbr'}))
            [out{:}] = reverse_bishop(P,V(revb),toli);
        elseif isfield(P,'Brev')
            [out{:}] = reverse_pvsyst(P,V(revb),toli,MaxIter);
        elseif all(isfield(P,{'Vbr','Vbi'}))
            [out{:}] = reverse_alonso(P,V(revb));
        else
            [out{:}] = forwardODM(P,V(revb),toli,MaxIter);
            warning('odm:revmodel','No reverse-bias parameters!');
        end
        I(revb) = out{1};
        if nargout > 1
            dIdV(revb) = out{2};
            d2IdV2(revb) = out{3};
        end
    end
end

function [i,di,d2i] = forwardODM(P,v,toli,MaxIter)
% Solve the 7-parameter diode model
%
%       i = Iph-Io·(exp(vj/nVth)-1)-vj/Rsh - Iph·di2mutau/(Vbi-vj)
%       vj = v + i·Rs

    % Get approx (Rs=0) recombination/reverse-voltage current
    Irec = P.Iph*P.di2mutau./(P.Vbi-v);
    
	% I = fiveparameterODM(V,P.Iph-Irec+Irev,P.Io,P.Rs,P.nVth,P.Rsh,toli,MaxIter);
    i = fiveparameterODM(v,P.Iph-Irec,P.Io,P.Rs,P.nVth,P.Rsh,toli,MaxIter);

    if (P.di2mutau > 0) && (P.Iph > 0)
        tricky = v > P.Vbi | Irec > exp(log(P.Io)+(v+i*P.Rs)/P.nVth);
        if any(tricky)
        % As V + IRs approaches Vbi, I -> -Inf, and V + IRs

            vv = v(tricky);
            K = 8;
            err = @(i) (P.Iph-exp(log(P.Io)+(vv+i*P.Rs)/P.nVth)+P.Io -(vv+i*P.Rs)/P.Rsh - ...
                       P.Iph*P.di2mutau./max(eps(0),P.Vbi-vv-P.Rs*i) - i)*K;

            i0 = P.Iph*(1-P.di2mutau/P.Vbi);
            i(tricky) = bisection(err,i0,-vv/P.Rs,@eps,@eps);
        end
    
        % Refine recombination/reverse-voltage currents (solve iteratively)
        for iter = 1:MaxIter+1
			oldI = i;
			Irec = P.Iph*P.di2mutau./max(eps(0),P.Vbi-v-P.Rs*i);
			i = fiveparameterODM(v,P.Iph-Irec,P.Io,P.Rs,P.nVth,P.Rsh,toli,MaxIter);
            if all(abs(i-oldI) <= toli(i)), break; end
        end
        if iter > MaxIter
            error('ODM could not converge to a solution'); 
        end
    end
       
	% Explicit 1st and 2nd derivatives
    if nargout > 1
        [di,d2i] = ODMderivatives(P,v + P.Rs*i);
    end
end

function [i,di,d2i] = reverse_bishop(P,v,toli)
% Bishop (1988) exponential model, solved by bisection, knowing that junction voltage must lie
% between v + i(vj = 0)·Rs and 0

    i0 = P.Iph*(1-P.di2mutau./P.Vbi);
    if P.Rs == 0
        i = v/P.Rsh*(1 + P.avalanche_a*(1-v/P.Vbr)^(-P.avalanche_m)); % Jsh
        i = P.Iph*(1-P.di2mutau./(P.Vbi-v)) - i;
        vj = v;
    else
        % jsh = @(vj) vj/P.Rsh.*(1 + P.avalanche_a*(1-vj/P.Vbr).^(-P.avalanche_m));
        err = @(vj) P.Iph*(1-P.di2mutau./(P.Vbi-vj)) - (vj-v)/P.Rs - vj/P.Rsh.*(1 + P.avalanche_a*(1-vj/P.Vbr).^(-P.avalanche_m));
        vj = bisection(err,max(P.Vbr,v+i0*P.Rs),0,@eps,toli);
        i = (vj - v)/P.Rs;
    end
    
    if nargout > 1
        K = 1./(P.Rsh*(P.Vbr - vj)).*P.avalanche_a.*(1-vj/P.Vbr).^(-P.avalanche_m);
        djA = K.*(P.Vbr + (P.avalanche_m-1)*vj);
        d2jA = K.*P.avalanche_m.*(2*P.Vbr + (P.avalanche_m-1)*vj)./(P.Vbr - vj);
        [di,d2i] = ODMderivatives(P,vj,djA,d2jA);
    end
end

function [i,di,d2i] = reverse_alonso(P,v)
% Alonso-García & Ruiz (2006) model, explicit in terms of V

    v(v < P.Vbr) = NaN;
    Be = 3.0;

    i0 = P.Iph*(1-P.di2mutau/P.Vbi);    % I @ Vj = 0
    m0 = ODMderivatives(P,0);           % dI/dV @ Vj = 0
    
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

function [i,di,d2i] = reverse_pvsyst(P,v,toli,MaxIter)
% PVsyst quadratic model: solution to I = Iph + ... + Brev(V+RsI)²

    if P.Rs == 0
        i = P.Iph*(1-P.di2mutau./(P.Vbi-v))-v/P.Rsh + P.Brev*v.^2;
    else
        i0 = P.Iph*(1-P.di2mutau/P.Vbi);    % I @ Vj = 0
        m0 = ODMderivatives(P,0);              % dI/dV @ Vj = 0
        
        % Start with the solution to i = i0 + m0(v+iRs) + Brev(v+iRs)²
        i = -v/P.Rs + 1/(2*P.Brev*P.Rs^2).*(1-sqrt(1-4*P.Brev*P.Rs*(1+P.Rs*m0)*(v + P.Rs*i0)));
        
        if P.di2mutau > 0
        % Fix deviations due to recombination current, if necessary  
            f = @(i) P.Iph*(1-P.di2mutau./(P.Vbi-v-i*P.Rs))-(v+i*P.Rs)/P.Rsh + P.Brev*(v+i*P.Rs).^2 - i;
            df = @(i) P.Rs*(-(P.Iph*P.di2mutau)./(P.Vbi-v-i*P.Rs).^2 -1/P.Rsh + 2*P.Brev*(v+i*P.Rs)) - 1;
            i = quickfzero(f,df,i,toli,MaxIter);
        end
    end
    
    if nargout > 1
        vj = v+i*P.Rs;
        djA = -2*P.Brev*vj;
        d2jA = -2*P.Brev;
        [di,d2i] = ODMderivatives(P,vj,djA,d2jA);
    end
end

function [dIdV,d2IdV2] = ODMderivatives(P,Vj,djA,d2jA)
% [dIdV,d2IdV2] = ODMDERIVATIVES(P,Vj,djA,d2jA) - Calculates first and second derivatives dI/dV
%   and d²I/dV² for the ODM curve at junction voltages Vj = V + Rs·I, optionally considering an 
%   avalanche-breakdown current Ia, or more specifically, its derivatives with respect to Vj:
%   djA = dIa/dVj, and d2jA = d²Ia/dVj² 

    if nargin < 3, djA = 0; end
    if nargin < 4, d2jA = 0; end
        
    djD = P.Io/P.nVth*exp(Vj/P.nVth);                % d(Idiode)/dVj
    djR = (P.Iph*P.di2mutau)./(P.Vbi-Vj).^2;         % d(Irec)/dVj
    djI = -(djD + djR + 1/P.Rsh + djA);              % dI/dVj
    dIdV = djI./(1 - P.Rs*djI);                      % dI/dV
    
    d2jI = -(2*djR./(P.Vbi-Vj)+djD/P.nVth + d2jA);   % d²I/dVj²
    d2IdV2 = (1+P.Rs*dIdV).^2./(1 - P.Rs*djI).*d2jI; % d²I/dV²
end
