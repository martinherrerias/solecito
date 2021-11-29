function setflag(id,flag,message,opt)
% SETFLAG(ID,FLAG,MESSAGE,OPT) - Change the status (enabled/color) and MESSAGE of the GUI step
%   ID, according to FLAG:
%
%       -3: Running (disabled)
%       -2: Disabled
%       -1: Error in last run (enabled)
%        0: Ready to run
%        1: Completed without errors
%        2: Completed with warnings - append 'Check warnings on log' after message
% 
%   indices are chosen so that: enabled = FLAG > -2, and completed = FLAG > 0
%
% MESSAGE: can be a string or cell array of strings (for multiple lines)
% OPT: 'app' Preserve existing text, add message as new line
%      'cat' Preserve existing text, concatenate message on last line
%       -N   (negative integer) Delete the last N lines and then append new message
%
% NOTE: when no message or opt are provided, default is setflag(id,flag,'','app'), i.e. preserve the
% existing message. Use setflag(id,flag,'') explicitly to clear all text.

    MAX_MSG_LEN = 500;
    global GUI
    prevmsg = get(GUI.(id).txt,'String');
    
    if nargin < 3, message = ''; opt = 'cat'; 
    elseif nargin < 4, opt = ''; end
    
    if ~iscell(message), message = cellstr(message); end
    if ~iscell(prevmsg), prevmsg = cellstr(prevmsg); end
    
    if isscalar(opt) && isnumeric(opt) && opt < 0 && mod(opt,1) == 0
       prevmsg = prevmsg(1:numel(prevmsg)+opt);
       opt = 'app';
    end
    
    switch opt
        case 'app', message = [prevmsg(:);message(:)];     
        case 'cat'
            prevmsg{end} = [prevmsg{end} message{1}];  
            message = [prevmsg(:);message(2:end)];
        case '' % do nothing
        otherwise
            error('setflag:opt','Unrecognized option: %s',opt)
    end
    
    warn = 'Check warnings on log';
    if flag == 2 && ~any(strcmp(message,warn)), message = [message;{warn}]; end
        
    % Check that messages have a sensible length
    toolong = cellfun(@numel,message) > MAX_MSG_LEN;
    message(toolong) = cellfun(@(s) [s(1:MAX_MSG_LEN - 3) '...'],message(toolong),'unif',0);
    
    set(GUI.(id).txt,'String',message);
    GUI.(id).flag = flag;
    
    switch flag
        case -3 % Running
            set(GUI.(id).txt,'ForegroundColor',[0 0 0.8]);
            set(GUI.(id).btn,'Enable','off');
        case -2 % Disabled
            set(GUI.(id).txt,'ForegroundColor',[0.6 0.6 0.6]);
            set(GUI.(id).btn,'Enable','off');
        case -1 % Error
            set(GUI.(id).txt,'ForegroundColor',[0.7 0 0]);
            set(GUI.(id).btn,'Enable','on');
        case 0 % Ready
            set(GUI.(id).txt,'ForegroundColor',[0 0 0]);
            set(GUI.(id).btn,'Enable','on');
        case 1 % Complete
            set(GUI.(id).txt,'ForegroundColor',[0 0.7 0]);
            set(GUI.(id).btn,'Enable','on');
        case 2 % Complete with warnings
            set(GUI.(id).txt,'ForegroundColor',[0.9 0.5 0]);
            set(GUI.(id).btn,'Enable','on');
        otherwise
            error('setflag:case','Say whaat?');
    end
    
    drawnow
end
