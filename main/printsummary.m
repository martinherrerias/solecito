function printsummary(simsteps,varargin)
% PRINTSUMMARY(SIMSTEP) - Print (updated) basic information about the step(s) SIMSTEP.
% The result should be a somewhat 'cleaner', yet similar output as that printed when executing
% a given step.

    global GUI
    
    nothing2print = @(x) any(GUI.(x).flag == [-2,-1,0]); % Disabled/Error/Ready-to-run
        
    % Parse
    allsteps = setdiff(fieldnames(GUI),{'menu'},'stable');
    if nargin < 1 || isempty(simsteps), simsteps = allsteps; end
    if ~iscell(simsteps),simsteps = {simsteps}; end
    assert(iscellstr(simsteps),'SIMSTEP must be a string, or cell-array of strings'); %#ok<ISCLSTR>
    unknown = setdiff(simsteps,allsteps);
    assert(isempty(unknown), shortliststr(unknown,'Unknown step','colon',':','quotes',''''));
 
    % Unless explicitly requested, stick to enabled & completed steps
    %flags = cellfun(@(s) GUI.(s).flag,simsteps);
    if nargin < 1
        simsteps(cellfun(nothing2print,simsteps)) = [];
        %flags = cellfun(@(s) GUI.(s).flag,simsteps);
    end
    
    for j = 1:numel(simsteps)
    % Calls to specific subfunctions
    
        % if nothing2print(simsteps{j}), printgeneric(simsteps{j}); continue; end
        
        switch simsteps{j}
            case 'setup', setupsummary(varargin{:});
            % case 'phystrck'
            % case 'physmod'
            % case 'layout'
            % case 'terrain'
            % case 'arrdef'
            % case 'reduce'
            % case 'meteo'
            % case 'shading'
            case 'irrtrans', irrtranssummary(varargin{:});
            % case 'models'
            case 'solver', solversummary(varargin{:});
                
            % PROVISIONAL:
            otherwise, printgeneric(simsteps{j});
        end
    end
end

function printgeneric(thestep)
% Generic Summary-Print (just copy GUI text)
    global GUI
    name = get(GUI.(thestep).btn,'String');
    msg = get(GUI.(thestep).txt,'String');
    if ~iscell(msg), msg = {msg}; end
    fprintf('%s (flag = %d):\n',name,GUI.(thestep).flag);
    cellfun(@(s) fprintf('\t%s\n',s),msg);
end

function setupsummary()
% Display (filtered), extended SimOptions structure
% FUTURE: if not running from the first time, display changes only

    global SimOptions;
    fprintf('\n');
    displayoptions(SimOptions);
    fprintf('\b');
end

function irrtranssummary(varargin)
% Print plant-averaged POA irradiance and shading results

    [opt,varargin] = getflagoptions(varargin,{'-verbose'});
    opt = getpairedoptions(varargin,opt,'restchk');

    EB = evalin('base','EB');
    reportmean = @(tag,v) fprintf(['\t' tag ': %0.1f kW/m²\n'],mean(v,'omitnan'));
    reportloss = @(tag,v,t) fprintf(['\t' tag ': %0.2f%%\n'],mean(~isnan(t).*v,'omitnan')/...
                                                             mean(~isnan(v).*t,'omitnan')*100);

    fprintf('\nRadiation balance (all-sample averages):\n\n');
    
    reportmean('Global Horizontal Irradiance',EB.GHI);
    if opt.verbose
        reportloss('\tFar-horizon & projection bias correction',EB.bias,EB.GHI);
        fprintf('\n');
    end
    
    Gpoa = EB.Bpoa0+EB.Dpoa0;

    reportloss('Tilt-gain',Gpoa - EB.GHI,EB.GHI);
    if opt.verbose
        reportmean('\tNon-shaded beam POA irradiance',EB.Bpoa0);
        reportmean('\tNon-shaded diffuse POA irradiance',EB.Dpoa0);
        reportmean('\tNon-shaded global POA irradiance',Gpoa);
        fprintf('\n');
    end
    
    reportloss('IAM factor (global)',EB.IAMd + EB.IAMb,Gpoa);
    if opt.verbose
        reportloss('\tIAM Loss (diffuse)',EB.IAMd,Gpoa);
        reportloss('\tIAM Loss (beam)',EB.IAMb,Gpoa);
        reportmean('\tIAM-adjusted global POA irradiance',Gpoa - EB.IAMd - EB.IAMb);
        fprintf('\n');
    end
    
    Gpoa = Gpoa - EB.IAMd - EB.IAMb;
    
    reportloss('Soiling Loss',EB.soiling,Gpoa);
    if opt.verbose, fprintf('\n'); end

    reportloss('Spectral Loss',EB.spectral,Gpoa);
    if opt.verbose, fprintf('\n'); end
    
    Gpoa = Gpoa - EB.soiling - EB.spectral;
    
    reportloss('Total shading Loss (irradiance)',EB.ShLb + EB.ShLd,Gpoa);
    if opt.verbose
        reportloss('\tDiffuse',EB.ShLd,Gpoa);
        reportloss('\tDirect',EB.ShLb,Gpoa);
        reportloss('\tDiffuse module mismatch (Bucciarelli)',EB.bucciarelli,Gpoa); 
        % reportmean('\tShaded global POA irradiance',Gpoa - EB.ShLd - EB.ShLb);
        fprintf('\t\t-----\n');
        fprintf('\t\tEstimated effect (with mismatch)\n');
        reportloss('\ta) Martínez-Moreno et al. 2010, central',EB.ShLaE,Gpoa);
        reportloss('\tb) Martínez-Moreno et al. 2010, string/mounts',EB.ShLaU,Gpoa);
        reportloss('\tc) Absolute worst case',EB.ShLaW,Gpoa);
        fprintf('\n');
    end

    reportmean('Global Effective Irradiance',EB.Ge);
    if opt.verbose
        reportmean('\tEffective Dpoa',EB.De);
        shf = nanmin(1 - EB.ShLb./EB.Bpoa0,1);
        reportmean('\tEffective Bpoa',EB.Be.*shf);
    end
    fprintf('\n');
end

function solversummary(varargin)
% Print plant-averaged electrical simulation results

    [opt,varargin] = getflagoptions(varargin,'-verbose');
    opt = getpairedoptions(varargin,opt,'restchk');

    % dt = evalin('base','Time.DoY');
    % dt = mode(round(diff(dt)*24*3600)/60); % minutes
                                                         
    reportsum = @(tag,v,u) fprintf(['\t' tag ': %0.2f %s\n'], mean(v,'all','omitnan'), u);
    reportloss = @(tag,v,t) fprintf(['\t' tag ': %0.2f%%\n'], mean(~isnan(t).*v,'all','omitnan')/...
                                                              mean(~isnan(v).*t,'all','omitnan')*100);
    UNIT = 'kW';
    MULT = 1e-3;
                                                          
    EB = evalin('base','EB');
    Mod = evalin('base','ModIVint');
    % POA = evalin('base','POA');
    ArrDef = evalin('base','ArrayDef');
    Trck = evalin('base','Trackers');
    MeteoData = evalin('base','MeteoData');
    celltemp = evalin('base','celltemp');
    
    Nu = numel(Trck.analysedtrackers);
    Nm = ArrDef.psize(2);
    
    allmodules = ArrDef.psubs(:);  % [mount-idx,module] indices for all connected modules
    [~,allmodules(:,1)] = ismember(allmodules(:,1),Trck.analysedtrackers);  % [used-mount-idx,..]
    allmodules = sub2ind([Nu,Nm],allmodules(:,1),allmodules(:,2));
    N = numel(allmodules);
    
    Pmp0 = Mod.Imp0*Mod.Vmp0;
    
    Ge = EB.Ge;
    Tc = celltemp(Ge,MeteoData.Ta,MeteoData.vw);
    %Ge = (1-POA.ShF).*POA.Bpoa0 + POA.Dpoa;
    % [~,ta,vw] = compatiblesize(Ge,MeteoData.Ta,MeteoData.vw);
    % Tc = celltemp(Ge,ta,vw);
    % clear ta vw
    % Ge = Ge(:,allmodules);
    % Tc = Tc(:,allmodules);
    
    fprintf('\nProduction balance (all-sample averages):\n\n');
        
    reportsum('Installed power',Pmp0*N*MULT,[UNIT,'p']);
    reportloss('Module efficiency at STC',Pmp0/(Mod.G0*Mod.area),1);
    
    reportsum('Nominal mean power @ STC efficiency',Ge/Mod.G0*Pmp0*N*MULT,UNIT);
    P0 = Ge/Mod.G0*Pmp0;
    Px = Mod.getMPP(Mod.G0,Tc,'MPPinterp').*Ge/Mod.G0;
    reportloss('\tVirtual loss due to temperature',P0 - Px,P0);
    
    P0 = Px;
    Px = Mod.getMPP(Ge,Mod.Tc0,'MPPinterp');
    reportloss('\tVirtual loss due to irradiance level',P0 - Px,P0);
    reportsum('Virtual power for non-shaded array',EB.E0*MULT,UNIT);
    fprintf('\n');

    reportloss('Shading & External-Mismatch Losses',EB.ShL,EB.E0);
    if opt.verbose
        reportloss('\tLinear (E)',EB.E0-EB.Elsh,EB.E0);
        reportloss('\tWorst-case (E)',EB.E0-EB.Ewcs,EB.E0);
        reportsum('Virtual mean power after shading',EB.Esh*MULT,UNIT);
    end
    
    fprintf('\n');
    reportloss('Module quality loss',0,1);
    reportloss('Intrinsic mismatch loss',0,1);
    reportloss('Ohmic wiring losses',0,1);
    reportloss('MPPT voltage limitation',EB.Esh-EB.Edc,EB.Esh);
    fprintf('\n');
    
    reportsum('Mean DC power at inverter(s) input',EB.Edc*MULT,UNIT);
    
    reportloss('Inverter Losses',EB.Edc-EB.Eac,EB.Edc);
    if opt.verbose
        reportloss('\tEfficiency',EB.inv.eff,EB.Edc);
        reportloss('\tPower threshold (lower)',EB.inv.clip.d,EB.Edc);
        reportloss('\tPower threshold (upper)', EB.inv.clip.p + EB.inv.clip.t,EB.Edc);
        reportloss('\tClipping (voltage)',EB.inv.clip.v,EB.Edc);
        reportloss('\tClipping (current)',EB.inv.clip.i,EB.Edc);
        reportloss('\tSelf consumption',EB.inv.self,EB.Edc);
    end

    fprintf('\n');
    reportsum('Mean converted AC Power',EB.Eac*MULT,strrep(UNIT,'W','VA'));
end

function displayoptions(options)
% Quick & Dirty function to display nested structure options

% TODO: (partly) replace by DISPNESTED?

    unimportant = {'DelimiterPriority'}; % SimOptions not to be included in the log file.
    
    fields = fieldnames(options);
    n = numel(fields);
    substruct = false(n,1);
    trash = false(n,1);
    
    for j = 1:n
        substruct(j) = isstruct(options.(fields{j}));
        if any(strcmp(fields{j},unimportant)), trash(j) = true; end
    end
    
    % Remove unimportant fields
    simpleopt = rmfield(options,fields(trash | substruct)); %#ok<NASGU>
    simpleopt = evalc('disp(simpleopt)'); 
    nospecialchars = @(x) regexprep(x,'([%\\])','$1$1'); % replace fprintf special characters
    fprintf(nospecialchars(simpleopt));
    
    sidx = find(substruct & ~trash);
    for j = 1:numel(sidx)
        s = options.(fields{sidx(j)});
        if isscalar(s)
            fprintf('\t%s:\n',fields{sidx(j)});
            displayoptions(s);
        else
            for k = 1:numel(s)
                fprintf('\t%s(%d):\n',fields{sidx(j)},k);
                displayoptions(s(k));
            end
        end
    end
end
