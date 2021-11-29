function POA = applyUF(POA,MeteoData,UFpars)

    if nargin < 3
        [FileName,PathName] = uigetfile('*.samlib','Utilization-Factor Parameter File');
        UFpars = pvlmod_SAMLibraryReader([PathName FileName]);
    end

    [Nt,Nat] = size(POA.Gpoa);
    AM = repmat(MeteoData.AMr,1,Nat);
    Ta = repmat(MeteoData.Ta,1,Nat);
    UF = zeros(Nt,Nat);
    UF(:) = UtilizationFactor(UFpars,AM(:),Ta(:),POA.Gpoa(:)/1000);

    POA.Gpoa_UF = POA.Gpoa;
    POA.Dpoa_UF = POA.Dpoa;
    
    POA.UF = reshape(UF,Nt,Nat);
    POA.Gpoa = POA.Gpoa.*UF;
    POA.Dpoa = POA.Dpoa.*UF;
end
