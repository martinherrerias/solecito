function [val, fromdefaults] = getSimOption(fieldnames,setifdefault,varargin)
% [VAL, FROMDEFAULTS] = GETSIMOPTION(FIELDNAME,SETIFDEFAULT)
%   Attempts to retrieve the SIMOPTIONS.FIELDNAME value from the global SIMOPTIONS structure.
%   If it fails, it tries to get the same value from DEFAULTOPTIONS().
%
%   FIELDNAME: strings of the form A.B.C... will look for subfields of structure elements.
%   SETIFDEFAULT: if true (default is false) and FIELDNAME is only found in DEFAULTOPTIONS, 
%       a matching field will be copied into the SimOptions structure. It's highly recommended to 
%       use COMPLETEOPTIONS instead, once and for all.
%   VAL: the value of the requested field
%   FROMDEFAULTS: is true if VAL was not set explicitly (i.e. it is the corresponding value from
%       the global OPTIONFLAGS structure).
%
% See also: GUISETUP, COMPLETEOPTIONS, DEFAULTOPTIONS

    if nargin < 2, setifdefault = false; end
    global SimOptions;
    global OptionFlags;
    
    persistent names;
    if isempty(names), names = {}; end
    
    if isempty(SimOptions)
        [SimOptions,OptionFlags] = completeoptions(); % contains call to DEFAULTOPTIONS
        warning('getSimOption:alldefaults','global SimOptions structure not found, setting to Defaults!');
    end
    if nargin < 1 || isempty(fieldnames), val = []; fromdefaults = false(0); return; end
    
    if iscell(fieldnames), S = struct(); D = struct(); end
    fields = cellstr(fieldnames);
    
    for j = 1:numel(fields)
        ff = strsplit(fields{j},'.');
        try
            val = getfield(SimOptions,{1},ff{:});
        catch
            [SimOptions,OptionFlags] = completeoptions();
            try
                val = getfield(SimOptions,{1},ff{:});
            catch
                try
                    match = parselist(fields{j},names,varargin{:});
                    assert(~isempty(match));
                catch
                    [explicit,implicit] = nestedfieldnames(SimOptions);
                    names = [explicit;implicit];
                    try
                        match = parselist(fields{j},names,varargin{:});
                        assert(~isempty(match));
                    catch
                        error('getSimOption:getdefault','SimOptions/Defaults.%s not found',fields{j});
                    end
                end
                if iscell(match)
                   error('getSimOption:ambiguous','Ambiguous SimOption match: %s',shortliststr(match));
                end
                ff = strsplit(match,'.');
                val = getfield(SimOptions,{1},ff{:});
                warning('getSimOption:ignorecase',...
                    'SimOptions.%s not found, using case-insensitive match %s',fields{j},match);
                fields{j} = match;
            end
        end
        
        if setifdefault
            SimOptions = setnestedfield(SimOptions,fields{j},val); 
            OptionFlags = setnestedfield(OptionFlags,fields{j},1); % set flag to 1 (default)
            warning('getSimOption:setdefault','Setting SimOptions.%s to Default',fields{j});
        end
        
        if nargout > 1, fromdefaults = getfield(OptionFlags,{1},ff{:}); end
        
        if j == 1 && ~iscell(fieldnames), return; end
        
        S = setfield(S,{1},ff{:},val);
        if nargout > 1, D = setfield(D,{1},ff{:},fromdefaults); end
    end
    
    val = S;
    fromdefaults = D;
    
%     if isnestedfield(SimOptions,fld)
%         val = getnestedfield(SimOptions,fld);
%     else
%         Defaults = DefaultOptions();
%         if isnestedfield(Defaults,fld)
%             val = getnestedfield(Defaults,fld);
%             warning('getSimOption:getdefault','SimOptions.%s not found, using Default',fld);
%             if setifdefault
%                 SimOptions = setnestedfield(SimOptions,fld,val); 
%                 OptionFlags = setnestedfield(OptionFlags,fld,1); % set flag to 1 (default)
%                 warning('getSimOption:setdefault','Setting SimOptions.%s to Default',fld);
%             end
%         else
%             error('getSimOption:getdefault','SimOptions/Defaults.%s not found',fld);
%         end
%     end
    
    % Check that there is a SimOption.xml file for MEX function calls
    % writeSimOptionXML();
end