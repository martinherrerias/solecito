function restorehandles()
% Restore UIControl object handles for buttons and text-boxes of GUI, either because they've been
% corrupted by mistake, or replaced by their text labels during saving.
%
% NOTE: RESTOREHANDLES assumes that Callback function names are @GUI<step> and that text-box- 
%   button pairs are the only UIControls in GUIfigure('GUImain').
%
% See also: GUIFIGURE

    global GUI

    h = get(GUIfigure('GUImain'),'children');
    h(arrayfun(@(x) ~isa(x,'matlab.ui.control.UIControl'),h))=[];
    txth = arrayfun(@(x) strcmpi(x.Style,'text'),h);
    h = [h(~txth),h(txth)];
    steps = arrayfun(@(x) regexprep(func2str(x.Callback),'^GUI(.*)$','$1'),h(:,1),'unif',0);

    for j = 1:numel(steps)
        txt = GUI.(steps{j}).txt;
        if ischar(txt) || iscellstr(txt) || isstring(txt)
            GUI.(steps{j}).txt = h(j,2);
            set(GUI.(steps{j}).txt,'String',txt);
        elseif ~ishandle(txt)
            GUI.(steps{j}).txt = h(j,2);
        end
        GUI.(steps{j}).btn = h(j,1);
    end
end