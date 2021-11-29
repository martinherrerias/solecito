function GUIarrdef(~,~)
% Read array definition *.arrdef file
% Requires existing layout structure Trackers (including Trackers.geom)

    global GUI;

    setflag('arrdef',-3,'Runing GUIarrdef...');
    finalwillandtestament = onCleanup(@updateflags);
    fprintf('\nRunning GUIarrdef...\n');
    
    Trck = evalin('base','Trackers');

    % Read tracker layout
    setflag('arrdef',-3,'Reading Array Definition...');
    fprintf('Reading Array Definition...\n');
    [arrdef,edims,pdims,StrIdx] = readarraydefinition(pickfile('*.arrdef'));

    % Check array definition
    physidx = {Trck.pidx,Trck.tidx};
    [arrdef,edims,pdims] = checkarraydefinition(arrdef,edims,pdims,physidx);

    % PROVISIONAL: ArrayDefinition only takes 5 columns (m,s,i,p,t), dump boxes to stridx
    % FUTURE: support for other coordinate types: MPPT, trafo, AC-busbar, block, etc.
    [arrdef,ArrIdx] = striparraydefinition(arrdef,edims,pdims,Trck.tidx,Trck.pidx);

    if ~isempty(StrIdx)
        ArrIdx = mergeidx(struct('is',StrIdx),ArrIdx,true);
    end
    
    Ntr = size(Trck.centers,2);
    Np = numel(Trck.geom.elements); % Number of Modules per Tracker
    
    psz = max(arrdef(:,[5,4]),[],1);
    if any(psz < [Ntr,Np])
        warning('GUIarrdef:pdimcheck','Some mounts/positions are not included in array definition');
    elseif any(psz > [Ntr,Np])
        error('GUIarrdef:pdimcheck', 'Array definition outside physical layout boundaries');
    end
    ArrayDef = NDmap(arrdef(:,[3,2,1]),arrdef(:,[5,4]),[],[Ntr,Np]);

    flagmsg{1} = sprintf('Complete Array Def.: %d mppts· %d strings · %d modules.',ArrayDef.esize);
    fprintf('\t%s\n',flagmsg{1});

    % TODO: is this necessary? isn't it already verified by checkarraydefinition?
        % Get No. of modules-per-string: for string k, in MPPT l, Nm = mps(k,l) 
        mps = shiftdim(sum(ArrayDef.PIDX>0,1),1); 
        % MPPTs can have different number of strings, i.e. mps(k,l) = 0, 
        % but not strings with different number of modules, mps(k,l) = {mps(q,l) or 0}
        problematic = any(mps > 0 & mps ~= repmat(mode(mps,1),ArrayDef.esize(2),1),1);
        if any(problematic)
            error(['You really want to check the Array Definition,',...
                    'there are strings with different numbers of modules in MPPT(s): {',...
                    repmat('%d ',1,nnz(problematic)) '}'],find(problematic));
        end

    % Get mount indices for all connected mounts
    Trck.analysedtrackers = unique(ArrayDef.psubs(:)*[1;0],'sorted');

    if ~runningfromUI()
        fprintf('\tNot-Interactive mode: not plotting\n');
    else
        if prod(ArrayDef.psize) < 10000
            plotting = true;
        else
            switch questdlg(['Plotting the complete array-definition can take a while, ',...
                             'sure you want to wait?'],'ArrayDef plotting')
                case 'Yes', plotting = true;
                otherwise, plotting = false;
            end
        end
        if ~plotting
            fprintf('\tArrayDef too large to plot...\n');
            flagmsg{2} = 'ArrayDef too large to plot';
            setflag('arrdef',-3,flagmsg);
        else
            fprintf('\tPlotting Array Definition...\n');
            setflag('arrdef',-3,{flagmsg{1},'Plotting Array Definition...'});
            fh = plotArrayDef(ArrayDef,Trck,ArrIdx);
            fprintf('\b\b\b\b: complete\n');
            if getSimOption('exportplots'), exportfigure('arraydef',fh); end
        end
    end

    assignin('base','ArrayDef',ArrayDef);
    assignin('base','ArrIdx',ArrIdx);
    assignin('base','Trackers',Trck);
    setflag('arrdef',1,flagmsg);

    % If Electrical-Models exist, check their compatibility with new ArrayDef
    if GUI.models.flag > 0
       if ~checkdesign(), setflag('models',2); else, setflag('models',1); end
    end
    
    % Run GUIreduce automatically on non-UI mode, or with non-deault opt.analysedmppts
    if ~strcmpi(getSimOption('analysedmppts'),'all') || ...
            ~runningfromUI() && ~isequal(Trck.analysedtrackers,1:Ntr)
        
        GUIreduce();
        
    elseif ~isequal(Trck.analysedtrackers(:)',1:Ntr)
    % Otherwise just issue a warning with instructions
        warning(['Only %d of %d mounts have connected modules, you can restrict any ',...
                'further analysis to those mounts by using menu Project > Reduce...'],...
                numel(Trck.analysedtrackers),Ntr);
    end
end
