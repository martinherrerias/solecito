function GUIirrtrans(~,~)
% GUIIRRTRANS - Wrapper for POAIRRADIANCE & EFFECTIVEIRRADIANCE, i.e. calculation of broad-band
%   effective plane-of-array irradiance at module level from ShadingResults and Meteo-Data (POA)
%   and plant-wide 'Energy Balance' summary statistics (EB)
%
%   THIS STEP CAN, AND PROBABLY SHOULD BE SKIPPED!: the resulting POA structure Nt·Nu·Nm can be
%   very large, and relatively cheap to recompute. i.e. it creates unnecessary large files and
%   when running GUISOLVER in batches (SimOption.runparallel = true), creates a big overhead when
%   splitting arguments. GUISOLVER can be called directly, and perform irradiance transposition
%   internally (at node/job level), keeping only EB and the final SolRes.

%     global GUI
    
    setflag('irrtrans',-3,'Runing GUIirrtrans...');
    finalwillandtestament = onCleanup(@updateflags);
    fprintf('\nRunning GUIirrtrans...\n');

    MeteoData = evalin('base','MeteoData');
    SunPos = evalin('base','SunPos');
    ShRes = evalin('base','ShRes');  
    Trck = evalin('base','Trackers');
    fHor = evalin('base','HorizonProfile.fHor');
    material = evalin('base','ModIVint.source.material');
    
%     Ntr = size(Trck.centers,2);
%     Nu = numel(Trck.analysedtrackers);
%     assert(any(ShRes.Nu == [Nu,Ntr]),'GUIirrtrans:Ntr','Inconsistent ShRes and Trck');

    fprintf('Evaluating Shading Results...\n');

    [POA,ShR] = poairradiance(MeteoData,SunPos,Trck,fHor,ShRes);
    [POA,EB] = effectiveirradiance(POA,MeteoData,Trck,material);

    flagmsg{3} = 'Irradiance-transposition complete, see details in log';

    assignin('base','POA',POA);
    assignin('base','EB',EB);
    assignin('base','ShRes_expanded',ShR);
    if getSimOption('resultsxls'), struct2csv({SunPos,EB},ShRes.Nt,'EnergyBalance.xls'); end

%     % If Array-Definition and Electrical-Models also exist, check them with new POA
%     if GUI.arrdef.flag > 0 && GUI.models.flag > 0
%        if ~checkdesign(), setflag('models',2); else, setflag('models',1); end
%     end
    
    printsummary('irrtrans','-verbose');
    setflag('irrtrans',1,flagmsg);
end