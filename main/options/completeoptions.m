function [options,fromdefaults] = completeoptions(options,defaults,varargin)
% OPTIONS = COMPLETEOPTIONS(OPT,DEFAULTS) - Complete the options structure OPTIONS with default 
%   values from structure DEFAULTS.
%
% [OPTIONS,FROMDEFAULTS] = COMPLETEOPTIONS(OPTIONS) - Use global SIMOPTIONS as the source of
%   default options, after an equivalent COMPLETEOPTIONS(SIMOPTIONS,DEFAULTOPTIONS)
%
% [OPTIONS,FROMDEFAULTS] = COMPLETEOPTIONS() - returns the set of current defaults, i.e. the
%   result of COMPLETEOPTIONS(SIMOPTIONS,DEFAULTOPTIONS)
%
% FROMDEFAULTS has the same structure as OPTIONS, but contains boolean field values indicating if 
%   any option OPTIONS.(j) comes from DEFAULTS.(j)

%
% See also: COMPLETESTRUCT

    global SimOptions
    global OptionFlags
    if isempty(SimOptions), SimOptions = struct(); OptionFlags = struct();
    elseif isempty(OptionFlags), OptionFlags = struct(); end
    
    if nargin < 1 || isempty(options), options = struct(); end
    if nargin < 2 || isempty(defaults) 
    % Complete global SimOptions with defaults (e.g. resolve version compatibility issues)
    
        % Remove obsolete options
        [val,names] = nestedstruct2cell(SimOptions);
        obsolete = ~isnestedfield(DefaultOptions(),names);
        if any(obsolete)
            warning('Removed obsolete %s',shortliststr(names(obsolete),'option'));
            SimOptions = cell2nestedstruct(val(~obsolete),names(~obsolete));
            
            [val,names] = nestedstruct2cell(OptionFlags);
            obsolete = ~isnestedfield(SimOptions,names);
            OptionFlags = cell2nestedstruct(val(~obsolete),names(~obsolete));
        end
    
        % Use defaults for any new options
        [SimOptions, isnew] = completestruct(SimOptions,DefaultOptions(),varargin{:});
        
        % Update global OptionFlags
        for f = nestedfieldnames(SimOptions)'
            if isnestedfield(OptionFlags,f{1})
                isdef = getnestedfield(OptionFlags,f{1});
            else
                isdef = getnestedfield(isnew,f{1});
            end
            OptionFlags = setnestedfield(OptionFlags,f{1},isdef);
        end
        
        defaults = SimOptions;
    end

    [options, fromdefaults] = completestruct(options,defaults,varargin{:});

    % get rid of incompatible option combinations
    if nargin < 2   % i.e. completing with defaults
        options.groundshading = options.groundshading && options.diffuseshading;
        options.anisotropic = options.anisotropic && options.diffuseshading;
    end
end