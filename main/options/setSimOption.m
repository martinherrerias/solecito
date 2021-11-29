function setSimOption(fieldname,val)
% SETSIMOPTION(FIELDNAME,VAL)
% Sets the value of SimOptions.fieldname = val, in global structure SimOptions
% Fieldaname strings of the form a.b.c... will create any missing nested substructures
% CAUTION: setting nested structure fields can create incomplete structures.
%
% See also: DEFAULTOPTIONS, GETSIMOPTION, GUISETUP

    global SimOptions;
    global OptionFlags;
    
    if isempty(SimOptions), getSimOption('version'); end  % force creation of defaults 

    if ~isnestedfield(SimOptions,fieldname)
        warning('setSimOption:newfield','%s is not already a member of SimOptions',fieldname);
    end
    SimOptions = setnestedfield(SimOptions,fieldname,val);
    OptionFlags = setnestedfield(OptionFlags,fieldname,0);  % Set flag to false (not default)
    
    % if ismember(fieldname,{'RelTol','minabstol','minreltol','NEPS','MaxIter'})
    % % Overwrite SimOption.xml file if any of its options are modified
    %     writeSimOptionXML(true);
    % end
end