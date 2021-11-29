function GUIlayout(~,~)
% Read mount coordinates (layout) definition from *.mounts file
% Requires pre-filled Trackers structure (phystrck)

    setflag('layout',-3,'\nRuning GUIlayout...');
    finalwillandtestament = onCleanup(@updateflags);
    fprintf('\nRunning GUIlayout...\n');
    
    % global GUI;

    Mounts = evalin('base','Trackers');

    % Read tracker layout
    setflag('layout',-3,'Reading layout...');
    fprintf('\nReading Layout...\n');
    Layout = importplantlayout(pickfile('*.mounts'));

    % FUTURE: Layout and Trackers should be separate structures, allowing the latter to be an
    % array of different Tracker types (models) and the former to specify their distribution (*)
    % PROVISIONAL: Merge them into Trackers structure
    Trackers = blendlayout(Layout, Mounts);

    Loc.latitude = Trackers.origin(2);
    Loc.longitude = Trackers.origin(1);

    Ntr = size(Trackers.centers,2);
    flagmsg{1} = sprintf('%s: %d %s-units',Trackers.name,Ntr,Trackers.type);
    fprintf('%s\n',flagmsg{1});

    setflag('layout',-3, [flagmsg,{'Checking mount clearance...'}]);

    allclear = all(checkmountclearance(Trackers) > 0);
    if ~allclear && getSimOption('groundshading')
        warning('checkgroundclearance:groundshading','Ground shades will be turned off')
        setSimOption('groundshading', false);
    end
    
    % Calculate tracker masks
    masktolerance = getSimOption('masktolerance');
    fprintf('\nCalculating tracker masks (%0.1e mask tolerance)...\n',masktolerance);
    Trackers.masks = gettrackermasks(Trackers,masktolerance);
    fprintf('%0.1f (avg.) neighbors/tracker\n',mean(sum(Trackers.masks,2)));

    % Plot
    setflag('layout',-3,[flagmsg,{'Plotting Layout...'}]);
    fh = plottrackerarray(Trackers);
    if getSimOption('exportplots'), exportfigure('layout',fh); end

    assignin('base','Trackers',Trackers);
    assignin('base','Location',Loc);
    
    setflag('layout',1,flagmsg);
end