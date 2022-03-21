function pvCplusplusSimOptions()
% PVCPLUSPLUSSIMOPTIONS() - make sure PVCPLUSPLUS_ROOT environmental variable (used by C++ code) 
%   is set to this function's working directory.
%
%   PROVISIONAL: keep PVCPLUSPLUS_ROOT/pvProj/Resources/SimOption.xml file updated with global
%   SimOptions
%   TODO: migrate to a common XML-based SimOptions file.

    RELPATH = 'pvProj/Resources/SimOption.xml';
    RELEVANT = {'RelTol','minabstol','minreltol','NEPS','MaxIter'};

    global SimOptions
    persistent lastOptions

    % make sure PVCPLUSPLUS_ROOT environment variable (used by C++ code) is set and works
    filepath = getenv("PVCPLUSPLUS_ROOT");
    if ~isempty(filepath), filepath = fullfile(filepath,RELPATH); end
    if ~isfile(filepath)
        filepath = fileparts(mfilename('fullpath'));
        setenv("PVCPLUSPLUS_ROOT",filepath)
        filepath = fullfile(filepath,RELPATH);
        assert(isfile(filepath),'Failed to find SimOption.xml file, dir. structure changed?');
    end
    
    % update SimOptions.xml file if something has changed
    if isempty(lastOptions), lastOptions = SimOptions; end
    
    if any(cellfun(@(f) SimOptions.(f) ~= lastOptions.(f),RELEVANT))
        lastOptions = SimOptions;
        writeSimOptionXML(SimOptions,filepath)
    end
end

function writeSimOptionXML(opt,file)

    txt = {'<?xml version="1.0" encoding="utf-8"?>';
    '<PrecTol>';
    sprintf('  <RelTol>%g</RelTol>',opt.RelTol);
    sprintf('  <minAbsTol>%g</minAbsTol>',opt.minabstol);
    sprintf('  <minRelTol>%g</minRelTol>',opt.minreltol);
    sprintf('  <NEPS>%d</NEPS>',opt.NEPS);
    sprintf('  <MaxIter>%d</MaxIter>',opt.MaxIter);
    sprintf('  <MaxRecur>%d</MaxRecur>',get(0,'recursionlimit'));
    '</PrecTol>'};

    fID = fopen(file,'w');
    cellfun(@(c) fprintf(fID,'%s\n',c),txt);
    fclose(fID);
end