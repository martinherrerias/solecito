function cleanlog(filename)
% CLEANLOG(FILE) Remove any intermediate <WAITBAR>(XXX.X%) ... lines from a log file.
% See also: OPTWAITBAR

    TAG = '<WAITBAR>'; % make sure it matches tag in OPTWAITBAR

    fID = fopen(filename,'r');
    lines =  textscan(fID,'%[^\n]');
    lines = lines{1};
    fclose(fID);

    matches = regexp(lines,['^' TAG '\(([\d.]+)%\)'],'tokens');
    toclear = ~cellfun(@isempty,matches);

    if any(toclear)
        matches = cat(1,matches{:});
        matches = str2double(cat(1,matches{:}));
        toclear(toclear) = matches > 0 & matches < 100;
    end

    if any(toclear)
        fID = fopen(filename,'w');
        fprintf(fID,'%s\n',lines{~toclear});
        fclose(fID);
    end
end