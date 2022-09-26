function [Pout,Loss] = ONDeval(Inv,Pin,Vin,Ta)
% [PAC,LOSS] = ONDEVAL(INV,PDC,VDC,TA) - Calculate the output of an inverter INV at input power 
%   PDC, input voltage VDC, and ambient temperature Ta, along with a detailed structure of losses.
%
% INPUT:
%   PDC - Input power
%   VDC - Input voltage (same size as PDC)
%   TA - Ambient temperature (same size as PDC)
%
%   INV.fPV - Inverter performance function Pac' = f(Pdc,Vdc), does not have to include clipping.
%   INV.MPPTHi/MPPTLow - MPPT operational voltage window.
%   INV.IdcMax - Absolute maximum input current
%   INV.Pdc0 - Absolute maximum input power (can be further limited by INV.TPLim)
%   INV.Ps0 - Minimum input power threshold
%   INV.Pnt/Paux - (optional) Night & Auxiliary losses, set to zero by default-.
%   INV.TPLim - (optional) Inverter temperature derate function f(Ta): Pac = min(Pac',f(Ta))
%       Default is @(x) INV.Pdc0, i.e. no temperature-dependent limit.
%
% OUTPUT:
%   POUT - size(PDC) array of output AC power values
%
%   LOSS.thresh - Power threshold losses, due to 0 < PDC < INV.Ps0
%   LOSS.self - Standby & night losses, due to INV.Pnt/Paux
%   LOSS.clip_i - Clipping losses due to max. input current
%   LOSS.clip_v - Clipping losses due to VDC outside MPPT window
%   LOSS.clip_p - Clipping losses due to PDC above 
%   LOSS.clip_t - Clipping losses due to temperature derate function
%   LOSS.eff - Losses due to inverter efficiency.
%
% See also: INVERTERPLOT, ONDREAD, PVLMOD_SAMLIBRARYREADER

    DEF.Pnt = 0;
    DEF.Paux = 0;
    DEF.TPLim = @(Ta)repmat(Inv.Pdc0,size(Ta));
    DEF.MPPTLow = 0;
    DEF.MPPTHi = Inf;
    DEF.IdcMax = Inf;

    Inv = completestruct(Inv,DEF);
    parsestruct(Inv,{'MPPTLow','MPPTHi','IdcMax','Ps0','Pdc0','Pnt','Paux'},'-n','-r','-p','-s');
    
    if isfield(Inv,'EffMax')
        validateattributes(Inv.EffMax,'numeric',{'scalar','real','finite','positive','<',1});
        if ~isfield(Inv,'fPV')
            Inv.fPV = @(Pdc,Vdc) max(0,Pdc - Inv.Ps0)*Inv.EffMax;
        end
    end

    parsestruct(Inv,{'TPLim','fPV'},'class',{'function_handle','griddedInterpolant'});
    
    if nargin < 3 || isempty(Vin) || all(isnan(Vin),'all')
       Vin =  0.5*(Inv.MPPTLow + Inv.MPPTHi);
       warning('Assuming Vdc = %0.1f',Vin);
    end

    [Pin,Vin,Ta] = compatiblesize(Pin,Vin,Ta);
    s = size(Pin);
    
    Pout = zeros(s);
    Loss.thresh = zeros(s);
    Loss.self = zeros(s);
    Loss.clip_i = zeros(s);
    Loss.clip_v = zeros(s);
    Loss.clip_p = zeros(s);
    Loss.clip_t = zeros(s);
    Loss.eff = zeros(s);
    
    dark = Pin < Inv.Ps0; 
    Loss.thresh(dark) = Pin(dark);
    Loss.self(dark) = Inv.Pnt + Inv.Paux;
    % Pout(dark) = -Inv.Pnt;
    % Pout(vout) = 0;
    
    % Assuming Solver already did its best effort to find an operational point in MPPT...
    vout = ~dark & ((Vin < Inv.MPPTLow) | (Vin > Inv.MPPTHi));
    Loss.clip_v(vout) = Pin(vout); 
    % Pout(vout) = 0;
    
    inrange = ~(dark | vout);
    
    pout = inrange & (Pin > Inv.Pdc0);
    Loss.clip_p(pout) = Pin(pout)-Inv.Pdc0;
    Pin(pout) = Inv.Pdc0;
    
    Iin = Pin./Vin;
    iout = inrange & Iin > Inv.IdcMax;
    Loss.clip_i(iout) = Pin(iout)-Inv.IdcMax*Vin(iout);
    Pin(iout) = Inv.IdcMax*Vin(iout);
        
    pout = inrange & (Pin > Inv.Pdc0);
    Loss.clip_p(pout) = Pin(pout)-Inv.Pdc0;
    
    Pout(inrange) = Inv.fPV(double(Pin(inrange)),double(Vin(inrange)));
    Loss.eff(inrange) = Pin(inrange) - Pout(inrange);
    
    Pmax = Inv.TPLim(Ta);
    tclip = Pout > Pmax;
    Loss.clip_t(tclip) = Pout(tclip)- Pmax(tclip);
    Pout(tclip) = Pmax(tclip);
end