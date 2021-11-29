function GUIreduce(analysedmppts)
% GUIREDUCE() - filter all existing structures, leaving only trackers connected to the given MPPTs
%
% PROVISIONAL: Operation defined by SimOptions.analysedmppts:
%   'all' - ignore analysedtrackers, mark everything for analysis
%   'trck' - set by analysedtrackers, issue warnings for incomplete MPPTs
%   [num. array] - ignore analysedtrackers, set by list of analysed MPPTs.
% Default (empty) is set to 'trck', or 'all' if isempty(analysedtrackers)

    global GUI;
    
    Trck = evalin('base','Trackers');
    ArrayDef = evalin('base','ArrayDef');
    ArrIdx = evalin('base','ArrIdx');
    HorProf = evalin('base','HorizonProfile');

    masktolerance = getSimOption('masktolerance');

    Ntr = size(Trck.centers,2);
    Ni = ArrayDef.esize(1);

    fprintf('Partial Analysis Set-up...\n');

    if nargin < 1 || ~any(analysedmppts)
    % PROVISIONAL: analysedmppts read from SimOptions file
        analysedmppts = getSimOption('analysedmppts');
    end

    if isempty(analysedmppts) || (ischar(analysedmppts) && strcmpi(analysedmppts,'trck'))
    % Get analysedmppts from analysedtrackers (info. only), issue warnings for incomplete MPPTs

        fprintf('Analysis bound by analysed-trackers\n');

        if ~isfield(Trck,'analysedtrackers') || ~any(Trck.analysedtrackers)
            warning('GUIreduce:Nat','No mounts are marked for analysis, setting them all on');
            Trck.analysedtrackers = 1:Ntr;
        end

        if ~isempty(analysedmppts)
        % Get a list of any MPPTs which are composed of analysedtrackers only
            mpptlist = unique(ArrayDef.esubs(Trck.analysedtrackers,:)*[1;0;0]); 
            allincluded = false(numel(mpptlist),1);
            for j = 1:numel(mpptlist)
                requiredmounts = unique(ArrayDef.psubs(mpptlist(j),:,:)*[1;0]);
                allincluded(j) = isempty(setdiff(requiredmounts,Trck.analysedtrackers));
            end
            analysedmppts = mpptlist(allincluded);
        end

        if isempty(analysedmppts)
            warning('GUIreduce:emptymppts',['SHADING ANALYSIS ONLY! No complete inverters in' ...
                    'analysed-mounts, shading analysis on the selected mounts is possible, ',...
                    'but it won''t work for an electrical simulation.']);
                
            setflag('solver',0,'No complete inverters in analysed-mounts');
        end

    elseif ischar(analysedmppts) && strcmpi(analysedmppts,'all')
    % Complete analysis: ignore analysedtrackers

        fprintf('Analysing all existing inverters\n');

        analysedmppts = (1:ArrayDef.esize(1))';
        %setSimOption('analysedmppts',analysedmppts);

        % Get mount indices for all connected mounts
        Trck.analysedtrackers = unique(ArrayDef.psubs(:)*[1;0]);

    else %if isnumeric(analysedmppts)
    % Explicit list of analysedmppts (internal numeric indices)

        fprintf('Analysing specific %s\n',shortliststr(analysedmppts,'MPPT'));

        % get internal (minimal) indices for analysedmppts
        n = numel(analysedmppts);
        if isfield(ArrIdx,'i')
            analysedmppts = ArrIdx.i.eidx(analysedmppts);
        else
            % ArrIdx.is.psize can change, if there are connection boxes:
            traildims = num2cell(ones(1,ArrIdx.is.Np - 1));
            analysedmppts = ArrIdx.is.esubs(analysedmppts,traildims{:})*[1;0];
        end

        assert(~isempty(analysedmppts),'GUIreduce:analysed',...
            'None of the provided analysedmppts could be found in the provided ArrayDef');
        if numel(analysedmppts) < n
            warning('GUIreduce:analysedeidx',...
            'At least one index in analysedmppts is not included in ArrayDef');
        end

        % Get mount indices that contain analysedmppts
        Trck.analysedtrackers = ArrayDef.psubs(analysedmppts,:,:);
        Trck.analysedtrackers = unique(Trck.analysedtrackers(:,1));
    end

    % Calculate tracker masks
    fprintf('\nCalculating tracker masks (%0.1e mask tolerance)...\n',masktolerance);
    Trck.masks = gettrackermasks(Trck,masktolerance);

    fprintf('%0.1f (avg.) neighbors/tracker\n',mean(sum(Trck.masks,2)));

    analysedtrck = Trck.analysedtrackers;
    Nat = numel(analysedtrck);
    Nam = numel(analysedmppts);

    % Remove trackers that are not analysed, nor are included in any mask
    [Trck,usedtrck,trckrelidx] = reducetrackers(Trck);
    if all(usedtrck)
       fprintf('All trackers are required, nothing to reduce\n\n');
       return;
    end
    
    Ntr = size(Trck.centers,2);
    Np = numel(Trck.pidx);
    
    fprintf('Reducing project to %d/%d mounts, %d/%d mppts\n',Nat,Ntr,Nam,Ni);

    % Update tracker layout plot, to show tracker masks, and possible reduction
    plottrackerarray(Trck);
    set(gcf,'Name','Reduced Mount Layout');

    % Reduce arrdef to analysedmppts
    arrdef = ArrayDef.deftable(analysedmppts,:,:);

    % Reduce to minimal indices. 
    % By passing trckrelidx.pidx, references to the new Trck numbering should be fixed
    [arrdef,arridx] = striparraydefinition(arrdef,'msi','pt',trckrelidx.pidx(:),Trck.pidx);

    % Merge relative references in arridx in the absolute index ArrIdx
    ArrIdx = mergeidx(ArrIdx,arridx);
    ArrDef = NDmap(arrdef(:,[3,2,1]),arrdef(:,[5,4]),[],[Ntr,Np]);

    if prod(ArrDef.psize) > 10000
        fprintf('\tReduced ArrayDef still too large to plot...\n');
        setflag('reduce',-3,flagmsg);
    else
        fprintf('\tUpdating Array Definition Plot...\n');
        plotArrayDef(ArrDef,Trck,ArrIdx);
        set(gcf,'Name','Reduced Array Definition');
    end
    
    assignin('base','Trackers',Trck);
    assignin('base','ArrayDef',ArrDef);
    assignin('base','ArrIdx',ArrIdx);
    
    % Reduce near-horizon-profiles
    if isstruct(HorProf) && ~isempty(HorProf.nHor)
        HorProf.nHor = HorProf.nHor(analysedtrck); 
        assignin('base','HorizonProfile',HorProf);
    end
    
    if GUI.shading.flag > 0
        % Reduce Shading-Results, if they were generated plant-wide
        ShRes = evalin('base','ShRes');
        ShRes = mountfilter(ShRes,analysedtrck);
        assignin('base','ShRes',ShRes);
    end

    if GUI.irrtrans.flag > 0 && evalin('base','exist(''POA'');')
        POA = evalin('base','POA');
        POA = filterstructure(POA,analysedtrck,'dim',2);
        assignin('base','POA',POA);
    end
    
    if GUI.solver.flag > 0
        
    end
end

function [trck,used,idx] = reducetrackers(trck)
% Reduce Trackers structure, removing not-used elements (not analysed, not included in masks)
  
    No = size(trck.centers,2); % Original size
    Na = numel(trck.analysedtrackers); % Analysed
    
    assert(all(size(trck.masks)==[Na,No]),'GUIreduce:reducetrck','Masks do not match analysedtrackers');

    % Boolean 1·No index 'used' includes analysed + any elements included in any mask
    used = any(trck.masks,1) | sparse(1,trck.analysedtrackers,true,1,No);
    if all(used), idx = []; return; end % Nothing to do
    
    % Remove rows/columns for ~used mounts
    trck.centers = trck.centers(:,used);
    if isfield(trck,'tilt'),trck.tilt = trck.tilt(used); end
    if isfield(trck,'azimuth'), trck.azimuth = trck.azimuth(used); end
    if isfield(trck,'slope'), trck.slope = trck.slope(used); end
   
    trck.masks = trck.masks(:,used);
    
    % re-make the index of analysedtrackers (note that this is a 'relative' index only, referred
    % to the sub-set of masked-trackers.
    idx = NDmap(find(used)');
    trck.analysedtrackers = idx.eidx(trck.analysedtrackers);
    
    % Update tracker index (refer trck numbers to absolute index)
    absidx = mergemap(NDmap(trck.tidx(:)),idx);
    trck.tidx = absidx.pidx(:);
end
