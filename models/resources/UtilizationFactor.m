function UF = UtilizationFactor(UFpars,AM,Ta,G)
% View Gerstmaier et al. 2011
% 'Validation of PVsyst Performance Model for the Concentrix CPV Technology'
%	AM - relative air mass
%	Ta - Ambient Temperature (ºC)
%	G - Irradiance (kW)

    if ~isstruct(UFpars)||~isnumeric(AM)||~isnumeric(Ta)||~isnumeric(AM)
        error('UtilizationFactor:InputParse','Something is wrong with input');
    end

    n = 0;
    if all(isfield(UFpars,{'wAM','AMthld','sAMlo','sAMhi'}))&&numel(AM)>0
        UFam = twolinefun(UFpars.AMthld,UFpars.sAMlo,UFpars.sAMhi,AM);
        n = numel(AM); wAM = UFpars.wAM;
    else, UFam = 0; wAM = 0;
    end
    if all(isfield(UFpars,{'wT','Tthld','sTlo','sThi'}))&&numel(Ta)>0
        UFt = twolinefun(UFpars.Tthld,UFpars.sTlo,UFpars.sThi,Ta);
        if n>0 && numel(Ta)~=n, numelerror; 
        else, n=numel(Ta); wT = UFpars.wT; end
    else, UFt = 0; wT = 0;
    end
    if all(isfield(UFpars,{'wE','Ethld','sElo','sEhi'}))&&numel(G)>0
       UFe = twolinefun(UFpars.Ethld,UFpars.sElo,UFpars.sEhi,G);
        if n>0 && numel(G)~=n, numelerror; 
        else, wE = UFpars.wE; end
    else, UFe = 0; wE = 0;
    end
    
    if wE+wT+wAM > 0
        UF = max(wAM*UFam + wT*UFt + wE*UFe,0);
    else
        UF = 1;
    end

    function y = twolinefun(thld,slo,shi,x)
        low = x<=thld;
        y = zeros(size(x));
        y(low) = 1+(x(low)-thld)*slo;
        y(~low) = 1+(x(~low)-thld)*shi;
    end
    function numelerror
       error('UtilizationFactor:InputParse','Vector Sizes must match');
    end
end

%%

% UFpars.wAM = 0.4; UFpars.AMthld = 1.7; UFpars.sAMlo = 0.1; UFpars.sAMhi = -0.1;
% UFpars.wT = 0.3; UFpars.Tthld = 25; UFpars.sTlo = 0.005; UFpars.sThi = 0;
% UFpars.wE = 0.3; UFpars.Ethld = 0.8; UFpars.sElo = 0.65; UFpars.sEhi = 0;
