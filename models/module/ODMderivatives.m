function [dIdV,d2IdV2] = ODMderivatives(P,Vj,djA,d2jA)
% [dIdV,d2IdV2] = ODMDERIVATIVES(P,Vj,djA,d2jA) - Calculates first and second derivatives dI/dV
%   and d²I/dV² for the ODM curve at junction voltages Vj = V + Rs·I, optionally considering an 
%   additional current term Ia (e.g. avalanche breakdown on reverse bias), or more specifically, 
%   its derivatives with respect to Vj: djA = dIa/dVj, and d2jA = d²Ia/dVj²
%
% See ONEDIODEMODEL for parameter description and application.

    if nargin < 3, djA = 0; end
    if nargin < 4, d2jA = 0; end
        
    djD = P.Io/P.nVth*exp(Vj/P.nVth);                % d(Idiode)/dVj
    djR = (P.Iph*P.di2mutau)./(P.Vbi-Vj).^2;         % d(Irec)/dVj
    djI = -(djD + djR + 1/P.Rsh + djA);              % dI/dVj
    dIdV = djI./(1 - P.Rs*djI);                      % dI/dV
    
    d2jI = -(2*djR./(P.Vbi-Vj)+djD/P.nVth + d2jA);   % d²I/dVj²
    d2IdV2 = (1+P.Rs*dIdV).^2./(1 - P.Rs*djI).*d2jI; % d²I/dV²
end
