function GUIsolver(~,~)
% GUISOLVER - Wrapper for PVARRAYSOLVER, core of the actual PV-system simulation (solution of DC 
%   electrical circuit); and ONDEVAL, the DC-AC conversion loss model.
%
%   Optionally, GUISOLVER can also replace GUIIRRTRANS, performing effective irradiance
%   transposition internally (at node/job level), keeping only the summary statistics EB and
%   discarding the bulky POA structure.

    global GUI;
    
    setflag('solver',-3,'Runing GUIsolver...');
    finalwillandtestament = onCleanup(@updateflags);
    fprintf('\nRunning GUIsolver...\n');
    
    ArrDef = evalin('base','ArrayDef');
    ModIVint = evalin('base','ModIVint');
    Diode = evalin('base','Diode');
    Inverter = evalin('base','Inverter');
    Trck = evalin('base','Trackers');

    celltemp = evalin('base','celltemp');
    MD = evalin('base','MD');
    
    [ArrDef,Trck] = connectedonly(ArrDef,Trck);
    
    if GUI.irrtrans.flag > 0
        try
            ShRes = evalin('base','ShRes_expanded');
            checkshadingresults(ShRes,ArrDef,Trck);
            POA = evalin('base','POA');
            EB = evalin('base','EB');
            assert(isequal(size(POA.Dpoa),[ShRes.Nt,ShRes.Nu,ShRes.Nm]),...
                'Inconsistent ShRes_expanded and POA');
        catch ERR
            setflag('irrtrans',-1,ERR.message);
            rethrow(ERR);
        end
    end
    if GUI.irrtrans.flag <= 0
        SunPos = evalin('base','SunPos');
        fHor = evalin('base','HorizonProfile.fHor'); 
        ShRes = evalin('base','ShRes');
        checkshadingresults(ShRes,ArrDef,Trck);
        if ~isequal(SunPos.Az,ShRes.az) || ~isequal(SunPos.El,ShRes.el)
            fprintf('\tInterpolating shading results...\n');
            ShRes = interpolate(ShRes,SunPos.Az,SunPos.El);
        end
    end
    
    parallel = getSimOption('runparallel');
    % M = getnworkers();
    % if M <= 1
    %     warning('Not enough memmory to span multiple workers');
    %     parallel = false; 
    % end
    
    [prjfolder,prjname] = fileparts(getSimOption('prjname'));
    resfilename = fullfile(prjfolder,['solres_' prjname '.mat']);
    [SolRes,flagmsg{1}] = checkforprevious(resfilename,ShRes.Nt,ArrDef);
    skipsolver = ~isempty(SolRes);

    if GUI.irrtrans.flag > 0
        if ~skipsolver
            [POA.Tush,POA.Tsh] = getcelltemperatures(POA,MD,celltemp);
            
            fprintf('\tRunning electrical solver...\n');
            setflag('solver',-3,[flagmsg(:),{'Running electrical solver...'}]);

            if parallel
            % Run DC solver on parallel batches
                [SolRes,cleaner] = runparallel(@pvArraySolver,...
                    {POA,ShRes,ArrDef,ModIVint,Diode,Inverter},1:2,[],...
                    'N',ShRes.Nt,'-simoptions','-backup');
            else
                % DEBUG Backdoor: Run single-threaded
                SolRes = pvArraySolver(POA,ShRes,ArrDef,ModIVint,Diode,Inverter);
            end
        end
        clear POA
    else 
        if skipsolver
            POA = poairradiance(MD,SunPos,Trck,fHor,ShRes);
            [~,EB] = effectiveirradiance(POA,MD,Trck,ModIVint.source.material);
            clear POA
        else
            fprintf('\tRunning combined transposition + solver...\n');
            setflag('irrtrans',-3,{'Running combined transposition + solver...'});
            setflag('solver',-3,{'Running combined transposition + solver...'});

            if parallel
            % Run POA irradiance + DC solver on parallel batches
                [SolRes,cleaner] = runparallel(@PVmain,...
                    {MD,SunPos,ShRes,Trck,fHor,ArrDef,ModIVint,celltemp,Diode,Inverter},...
                    1:3,[],'N',ShRes.Nt,'-simoptions','-backup');

                % PROVISIONAL: EB is just the result of EFFECTIVEIRRADIANCE, packed into SolRes
                % until runparallel takes nargout > 1     
                EB = SolRes.EB; SolRes = rmfield(SolRes,'EB');       
            else
                [SolRes,EB] = PVmain(MD,SunPos,ShRes,Trck,fHor,ArrDef,ModIVint,...
                    celltemp,Diode,Inverter);
            end
        end
    end
    save(resfilename,'SolRes'); % provisional?
    fprintf('\tSolver-Results saved to: %s\n',resfilename);

    timestr =@(x) [num2str(floor(x/3600)) ':' datestr(x/(24*3600),'MM:SS')];
    fprintf('\tSolver finished! (%s elapsed)\n',timestr(SolRes.simtimer.globaltime));

    % Apply DC-AC conversion losses
    [SolRes,L] = modelinverter(SolRes,Inverter,MD);
    EB = completestruct(EB,L);
    
    if getSimOption('resultsxls'), struct2csv(EB,ShRes.Nt,'SimResults.xls'); end

    assignin('base','SolRes',SolRes);
    assignin('base','EB',EB);

    if GUI.irrtrans.flag <= 0
        setflag('irrtrans',-1,'(completed within solver)'); 
        printsummary('irrtrans','-verbose');
    end
    setflag('solver',1,'Solver finished!');
    printsummary('solver','-verbose');
    
    fsave = get(GUI.menu.save,'Callback'); fsave();     % auto-save project
    
    if exist('cleaner','var'), cleaner(); end % if all went well, delete any partial results
end

function [ArrDef,Trck] = connectedonly(ArrDef,Trck)
    Nu = size(Trck.centers,2);  
    if isfield(Trck,'analysedtrackers') && ~isequal(Trck.analysedtrackers,1:Nu)
        
        [elist,plist] = map2lists(ArrDef);
        idx = NDmap(Trck.analysedtrackers);
        plist(:,1) = idx.eidx(plist(:,1));
        ArrDef = NDmap(elist,plist);
        
        Trck = filterstructure(Trck,Trck.analysedtrackers,Nu);
        Trck = filterstructure(Trck,Trck.analysedtrackers,Nu,'dim',2);
        Trck.analysedtrackers = 1:numel(Trck.analysedtrackers);
    end
end

function checkshadingresults(ShRes,ArrDef,Trck)

    assert(all([Trck.geom.dims(1),ArrDef.psize(2)] == ShRes.Nm),...
        'Inconsistent Nm in ShRes, ArrDef, and/or Trackers');
    assert(all([numel(Trck.analysedtrackers),ArrDef.psize(1)] == ShRes.Nu),...
        ['Inconsistent Nu in ShRes, ArrDef, and/or Trackers. Try reloading them using',...
         'GUIlayout, GUIarrdef, and GUIshading (in that order)']);
end

function varargout = PVmain(varargin)
% PVMAIN(MD,SunPos,ShRes,Trck,fHor,ArrDef,ModIV,celltemp,Diode,Inverter,[OPT,'name',val])
% PVMAIN('backup',S) - attempt to resume interrupted calculation from backup structure S,
%   where S = load(FILE), after a previous call with ..,'backup',FILE.

    global SimOptions
    SimOptions = completeoptions(); % fill with defaults if necessary
    options = SimOptions;
    
    options.backup = '';
    [options,varargin] = getpairedoptions(varargin,options);
   
    switch numel(varargin) 
    case 0
    % Resume from crash (already on pvArraySolver)
        assert(isstruct(options.backup) && isfield(options.backup,{'opt','POA'}),...
            'Failed to recover from backup')
        B = options.backup;
        options = rmfield(options,'backup');
        
        % New SimOptions can override backup options!
        options = completestruct(options,B.opt);
        
        fld = setdiff(fieldnames(B.opt),'backup');
        changed = ~cellfun(@(f) isequal(B.opt.(f),options.(f)),fld);
        if any(changed)
            warning('Attempting to resume with altered simulation %s',...
                shortliststr(fld(changed),'option','colon',':'));
            for j = find(changed), B.opt.(fld{j}) = options.(fld{j}); end
        end

        SolRes = pvArraySolver('backup',B.backup);
        [~,EB] = effectiveirradiance(B.POA,B.MD,B.Trck,B.ModIV.source.material);
        clear B
        
    case {10,11}
    % MD,SunPos,ShRes,Trck,fHor,ArrDef,ModIV,Diode,Inverter,[OPT])
        if numel(varargin) == 11
            options = completestruct(varargin{end},options); 
        end
        [MD,SunPos,ShRes,Trck,fHor,ArrDef,ModIV,celltemp,Diode,Inverter] = deal(varargin{1:10});
        
        fprintf('Evaluating Shading Results...\n');
        [POA,ShRes] = poairradiance(MD,SunPos,Trck,fHor,ShRes);
        [POA,EB] = effectiveirradiance(POA,MD,Trck,ModIV.source.material);

        [POA.Tush,POA.Tsh] = getcelltemperatures(POA,MD,celltemp);

        fprintf('\tRunning electrical solver...\n');
        SolRes = pvArraySolver(POA,ShRes,ArrDef,ModIV,Diode,Inverter,options);
    end
    
    if nargout > 1
        varargout = {SolRes,EB};
    else
    % runparallel requires single output, to stuff everything together 
        SolRes.EB = EB; 
        varargout = {SolRes};
    end
end

function [SolRes,EB] = modelinverter(SolRes,Inverter,MD)
    
    SolRes.Pac = zeros(size(SolRes.Vdc));
    Pdc = sum(SolRes.Pdc,3,'omitnan');
    Pdc(Pdc<0) = 0;
    SolRes.Vdc(SolRes.Vdc<0) = 0;
    [SolRes.Pac(:),SolRes.Loss] = ONDeval(Inverter,Pdc,SolRes.Vdc,MD.Ta);

    EB.E0 = sum(SolRes.Push,2:3,'omitnan');
    EB.Elsh = sum(SolRes.Plsh,2,'omitnan');
    EB.Ewcs = sum(SolRes.Pwcs,2,'omitnan');
    EB.Esh = sum(SolRes.Psh,2:3,'omitnan');
    EB.Edc = sum(SolRes.Pdc,2:3,'omitnan');
    EB.Eac = sum(SolRes.Pac,2,'omitnan');
    EB.ShL = EB.E0-EB.Esh;

    EB.inv.self = sum(SolRes.Loss.self,2,'omitnan');
    EB.inv.eff = sum(SolRes.Loss.eff,2,'omitnan');
    EB.inv.clip.d = sum(SolRes.Loss.thresh,2,'omitnan');
    EB.inv.clip.v = sum(SolRes.Loss.clip_v,2,'omitnan');
    EB.inv.clip.i = sum(SolRes.Loss.clip_i,2,'omitnan');
    EB.inv.clip.p = sum(SolRes.Loss.clip_p,2,'omitnan');
    EB.inv.clip.t = sum(SolRes.Loss.clip_t,2,'omitnan');

    EB.inv.clip.v = EB.inv.clip.v + sum(EB.Esh-EB.Edc,2,'omitnan');
    %EB.inv = sum(EB.Edc-EB.Eac)/sum(EB.Edc);
end

function [Tush,Tsh] = getcelltemperatures(POA,MD,celltemp)
% Get cell temperatures for shaded and not shaded elements (before applying UF!)

    fprintf('\tCalculating cell-temperatures...\n');
    Gush = POA.Dpoa + POA.Bpoa;
    [Gush,ta,vw] = compatiblesize(Gush,MD.Ta,MD.vw);
    Tush = celltemp(Gush,ta,vw);  
    Tsh = celltemp(POA.Dpoa,ta,vw);
end

function [SolRes,msg] = checkforprevious(resfilename,Nt,ArrDef)
% Look for existing results, to skip pvArraySolver when e.g. just inverter model changes
% TODO: should this even exist (once the code is stable)? AC modeling should be moved to a 
%   different simulation step, and in that case, under which real-use scenario would you want to 
%   reuse results?
%   If it stays, we need way stronger conformity checks, maybe turn solver-results into a class
%   with enough metadata, like SHADINGRESULTS.

    Ni = ArrDef.esize(1);                       % No. of inverters (MPPT-inputs)
    Ns = ArrDef.esize(2);                       % Max. No. of strings per MPPT 

    SolRes = [];
    if isfile(resfilename)
        switch optquestdlg('Solver results seem already available, do you want to use them?')
        case 'Yes'
            fprintf('\tChecking precalculated solver results...');
            setflag('solver',-3,'Checking precalculated solver results...');
            load(resfilename,'SolRes');
            fullsize = @(x) arrayfun(@(j) size(x,j),1:3);
            if ~all(fullsize(SolRes.Pdc)==[Nt,Ni,Ns]) % || some other conditions ...
                fprintf('sorry, they just do not seem right\n');
                SolRes = [];
                msg = sprintf('Pre-calculated results could not be used,');
            else
                fprintf('they look allright!\n');
                msg = 'Using precalculated solver results.';
            end
        case 'No'
            msg = 'Not using precalculated solver results,';
        otherwise
            error('GUIsolver:questdlg','You might want to try ''Yes'' or ''No'' next time')
        end
    else
        fprintf(['\tNo precalculated solver results found. Rename an old/modified ',...
                  'result-file as: <%s> if you want it to be found\n'],resfilename);
        msg = 'No precalculated solver results found,';
    end
end