function varargout = restoreSimOptions(backup)
% C = RESTORESIMOPTIONS() - create an internal copy of global SimOptions and OptionFlags, and
%   return an ONCLEANUP object C that, when destroyed, will set both variables back to their
%   original values. Use it to make provisional changes, and ensure they are reverted even if the
%   calling function crashes.
%
% EXAMPLE:
%   getSimOption('version')
%   tmp = restoreSimOptions();
%   setSimOption('version','foo');
%   getSimOption('version')
%   clear('tmp')
%   getSimOption('version')

    narginchk(0,1);
    global SimOptions;
    global OptionFlags;

    if nargin == 0 && nargout == 1
        backup = struct('SimOptions',SimOptions,'OptionFlags',OptionFlags);
        varargout{1} = onCleanup(@(x) restoreSimOptions(backup));
    elseif nargin == 1 && nargout == 0
        SimOptions = backup.SimOptions;
        OptionFlags = backup.OptionFlags;
    end
end