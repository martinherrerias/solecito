function GUIphystrck(~,~)
% Read physical mount definition from *.tpoly file
% Requires pre-filled Trackers structure (setup)

    global GUI

    setflag('phystrck',-3,'Runing GUIphystrck...');
    finalwillandtestament = onCleanup(@updateflags);
    fprintf('\nRunning GUIphystrck...\n');
            
    if GUI.physmod.flag > 0
        module = evalin('base','Trackers.module'); % salvage existing module
    end

    % Tracker-Area Definition {Trackers.geom}  
    fprintf('Reading physical tracker definition...\n');
    Trackers = importpolygonsfile(pickfile('*.tpoly','Pick tracker-polygons file'));

    fprintf('\tMount ID: %s\n',Trackers.name);

    [w,h] = rectangleproperties(Trackers.geom.border);

    Trackers.border = Trackers.geom.border;
    Nm = Trackers.geom.dims(1);

    type = typestring(Trackers.type);
    flagmsg = sprintf('%d modules, (%0.1f x %0.1f m) %s',Nm,w,h,type);
    fprintf('\t%s\n',flagmsg);

    assignin('base','Trackers',Trackers);

    figh = GUIfigure('trckareas','Physical Mount'); clf(figh);
    labeling = numel(Trackers.geom.elements)<100;
    plotArea(Trackers.geom,2,labeling);
    set(gca,'Color',[0.8;0.8;0.8]); axis equal; axis tight; title('Tracker'); xlabel('x [m]'); ylabel('y [m]');
    hold on
    plot(Trackers.analysedpoints(1,:),Trackers.analysedpoints(2,:),'kx','MarkerSize',10)

    if getSimOption('exportplots'), exportfigure('phystrck',figh); end

    setflag('phystrck',1,{'Physical Tracker Defined',flagmsg});
 
    if GUI.physmod.flag > 0     
        updateflags();
        GUIphysmod(module);
    end
end

function str = typestring(type)
    switch type
        case '0a' , str = 'fixed-tables';
        case '1aC', str = '1a-contour-trackers';
        case '1aF', str = 'fixed-tilt-1a-trackers';
        case '1aV', str = 'vertical-1a-trackers';
        case '2a' , str = '2-axis-trackers';
    end
end
