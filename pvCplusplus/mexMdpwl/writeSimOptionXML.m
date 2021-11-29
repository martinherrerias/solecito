
function writeSimOptionXML(overwrite)
% PROVISIONAL: check if there is a SimOption.xml file in the current directory, for calls
%   to MEXADDSERIES and MEXADDPARALLEL. If not (or if OVERWRITE == TRUE), write a new one
%   with the contents from global SimOptions.

    if nargin < 1 || isempty(overwrite), overwrite = false; end
    if ~overwrite && ~isempty(dir('SimOption.xml')), return; end

    global SimOptions

    txt = {'<?xml version="1.0" encoding="utf-8"?>';
    '<PrecTol>';
    sprintf('  <RelTol>%g</RelTol>',SimOptions.RelTol);
    sprintf('  <minAbsTol>%g</minAbsTol>',SimOptions.minabstol);
    sprintf('  <minRelTol>%g</minRelTol>',SimOptions.minreltol);
    sprintf('  <NEPS>%d</NEPS>',SimOptions.NEPS);
    sprintf('  <MaxIter>%d</MaxIter>',SimOptions.MaxIter);
    sprintf('  <MaxRecur>%d</MaxRecur>',get(0,'recursionlimit'));
    '</PrecTol>'};

    fID = fopen('SimOption.xml','w');
    cellfun(@(c) fprintf(fID,'%s\n',c),txt);
    fclose(fID);
end