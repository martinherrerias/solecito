function varargout = PV_simple(plant,data,mode,varargin)
% [opt,MBE,RMS,MAD] = PV_SIMPLE(PLANT,DATA,MODE,[PARAMS,MODULE,INVERTER,..])
% INPUT:
%   DATA - TIMETABLE with fields {Pdc,Pac,sunaz,sunel,Ta,vw,Tm}
%   MODEL - 'mattei' or 'odm+inv'

    narginchk(3,Inf);
    
    mode = parselist(mode,{'fit','run'}');
    validateattributes(plant,'struct',{'scalar'});
   
    opt.params = [];
    opt.Module = [];
    opt.Inverter = [];
    
    % Defaults
    opt.model = [];
    opt.pvModType = 'polysi';
    opt.albedo = 0.2;
    opt.CSR = 25;
    opt.IAM_maxloss = 0.0622;
    opt.mgap = 0.02;
    opt.BPR = 3;
    opt.strtypes = [];
    opt.SPI = [];
    opt.maxRPS = 1; % max rows per string
    
    opt.mismatch = 'martinezmoreno';
    opt.mismatch_opt = {};
    
    % fitting
    opt.knownGTI = false;
    opt.fitinverter = [];
    opt.fitSPI = [];
    opt.normalizeSPI = true;
    opt.fitmethod = 'lsq'; % lar, lgsq
    opt.TmWeight = 5e-3;
    opt.TmRow = [];
    
    % filters
    opt.min_GTI = 5;
    opt.min_sunel = 2;
    opt.max_eta = 0.3;
    opt.min_Pdc = 0.01;
    opt.min_Pac = 0.01;
    opt.clip_Pdc = 0.95;
    opt.max_Pac = 1.5;

    FLD.opt = {'params','Module','Inverter','model',...
               'pitch','tilt','azimuth','mdims','rows','mgap','BPR',...
               'PmpRef','PacMax','MPS','SPI','nSPI','strtypes','maxRPS',...
               'mismatch','mismatch_opt',...
               'pvModType','albedo','IAM_maxloss','CSR',...
               'fitmethod','knownGTI','fitinverter','fitSPI','normalizeSPI','TmWeight','TmRow',...
               'min_GTI','min_sunel','max_eta','min_Pdc','min_Pac','clip_Pdc','max_Pac',...
               'plot','verbose','dryrun'};

    % Parse PLANT, merge with OPT into a complete list of parameters
    plant = completestruct(plant,opt,'valid',@(x) ~isempty(x));
    [opt,isdef] = completestruct(plant,opt,'nested',false);
    unknown = setdiff(fieldnames(opt),FLD.opt);
    if ~isempty(unknown)
       warning(shortliststr(unknown,'Unknown field'));
       opt = rmfield(opt,unknown);
    end
    opt = orderfields(opt,intersect(FLD.opt,fieldnames(opt),'stable'));
    
    [wasdef,fld,implicit] = nestedstruct2cell(isdef);
    for j = 1:numel(implicit)
        s = fld(contains(fld,implicit{j}));
        fld(s) = [];
        fld(end+1) = implicit(j);
        wasdef(s);
    end
    
    [opt,~,isdef] = parseoptions(varargin,{'-plot','-verbose','-dryrun'},opt,'dealrest',3);
    for j = find(~cat(2,wasdef{:})), isdef.(fld{j}) = false; end
    
    if isempty(opt.model)
       if isempty(opt.Module), opt.model = 'mattei'; else, opt.model = 'odm+faiman'; end
    end 
    opt.model = parselist(opt.model,{'mattei','odm+faiman'}');
    opt.fitmethod = parselist(opt.fitmethod,{'lsq','lar','lgsq'});
    
    % Check mismatch options
    [~,foo] = mismatch(1,10,0.2,'model',opt.mismatch,opt.mismatch_opt{:});
    opt.mismatch = foo.model;
    
    if isempty(opt.params), opt.params = modeldefaults(opt.model); end
    
    % Parse DATA, remove unused variables and invalid rows
    FLD.vec = {{'sunaz','sunel','Ta','vw','DHI','BNI','ENI'},{'Tm','GTI','AMa','tpw'}}; % {req,opt}
    FLD.mat = {{'Pac','Pdc'},{'Vdc','Idc'}};
    switch mode
    case 'fit', FLD.mat = {{'Pac','Pdc'},{'Vdc','Idc'}};
    case 'run', FLD.mat = {{},{'Pac','Pdc','Vdc','Idc'}};
    otherwise, error('this should not happen');
    end

    if ~opt.dryrun || ~isempty(data)
        validateattributes(data,'timetable',{'nonempty'});
        foo = setdiff(data.Properties.VariableNames,[FLD.vec{:},FLD.mat{:}]);
        data(:,foo) = [];
        data = timetable2struct(data);

        if isfield(data,'Pac'), Ni = size(data.Pac,2);
        elseif isfield(data,'Pdc'), Ni = size(data.Pdc,2);
        else, Ni = NaN;
        end

        parsestruct(data,FLD.vec{1},'opt',FLD.vec{2},'real','vector');

        parsingpending = isnan(Ni);
        if ~parsingpending
            parsestruct(data,FLD.mat{1},'opt',FLD.mat{2},'real','2d','size',[NaN,Ni]);
        end

        useful = data.sunel > opt.min_sunel & isfinite(data.Ta) & isfinite(data.vw) & ...
            isfinite(data.BNI) & isfinite(data.DHI);

        if opt.knownGTI
            useful = useful & data.GTI > opt.min_GTI;
        elseif isfield(data,'GTI'), data = rmfield(data,'GTI');
        end
        data = filterstructure(data,useful);

        % Fill missing fields with []
        foo = setdiff([FLD.vec{:},FLD.mat{:}],fieldnames(data));
        for j = 1:numel(foo), data.(foo{j}) = zeros(numel(data.sunel),0); end
        
        if isempty(data.Tm), opt.TmWeight = 0; end
    else
        Ni = NaN;
        parsingpending = false;
    end
    
    % Check Inverter / Module
    opt.Module = opt.Module;
    if ~isempty(opt.Inverter)
        validateattributes(opt.Inverter,'struct',{'scalar'});
        ONDeval(opt.Inverter,0,0,0);
        if isdef.fitinverter, opt.fitinverter = false; end
    else
        opt.Inverter = struct.empty;
        if isdef.fitinverter, opt.fitinverter = strcmp(mode,'fit'); end
    end
    if ~isempty(opt.Module)
        if isstruct(opt.Module), opt.Module = ModuleInterpolant(opt.Module); end
        validateattributes(opt.Module,'ModuleInterpolant',{'scalar'});
    else
        assert(strcmp(opt.model,'mattei'),'odm+faiman requires provided ModuleInterpolant');
    end
    
    parsestruct(opt,{'pitch','tilt','azimuth','rows','PmpRef','PacMax','MPS'},...
                'opt',{'mgap','BPR'},'numeric','scalar','real');
            
    parsestruct(opt,{'pitch','PmpRef','PacMax'},'numeric','scalar','real','positive');
    parsestruct(opt,{'mdims'},'numeric','vector','real','positive','numel',2);
    parsestruct(opt,{'tilt'},'opt',{'mgap'},'numeric','scalar','real','nonnegative');
    parsestruct(opt,{'rows','MPS'},'opt',{'BPR'},'numeric','scalar','real','integer','positive');
    
    if isfield(opt,'nSPI')
        parsestruct(opt,{'nSPI'},'numeric','real','integer','positive','vector');
        if isnan(Ni), Ni = numel(opt.nSPI); end
        opt.nSPI = opt.nSPI(:)';
    end
    if isempty(opt.strtypes)
    % Create every combination of up to maxRPS rows
        opt.strtypes = cell(opt.maxRPS,1);
        for j = 1:opt.maxRPS
            opt.strtypes{j} = num2cell(nchoosek(1:opt.rows,j),2);
        end
        opt.strtypes = cat(1,opt.strtypes{:})';
    end
    if iscell(opt.strtypes)
        cellfun(@(x) validateattributes(x,{'numeric'},{'real','integer','positive','<=',opt.rows,'size',[1,NaN]}),opt.strtypes);
        jj = arrayfun(@repelem,1:numel(opt.strtypes),cellfun(@numel,opt.strtypes),'unif',0);
        opt.strtypes = full(sparse([jj{:}],[opt.strtypes{:}],true,numel(opt.strtypes),opt.rows));
    end
    validateattributes(opt.strtypes,{'numeric','logical'},{'2d','binary','size',[NaN,opt.rows]})
    
    u = size(opt.strtypes,1);
    if isempty(opt.SPI)
        assert(isfield(opt,'nSPI'),'nSPI and/or SPI must be provided');
        
        opt.SPI = repmat(ones(1,u)/u,Ni,1).*opt.nSPI'; % equal string-type "weights" for each inverter
        if isdef.fitSPI, opt.fitSPI = true; end
        isdef.SPI = true;
    else
        if size(opt.SPI,1) == 1 && Ni > 1, opt.SPI = repmat(opt.SPI,Ni,1); end
        
        if ~isfield(opt,'nSPI')
            validateattributes(opt.SPI,'numeric',{'real','integer','nonnegative','2d','size',[Ni,u]});
            opt.nSPI = sum(opt.SPI,2)';
            if isnan(Ni), Ni = numel(opt.nSPI); end
        else
            validateattributes(opt.SPI,'numeric',{'real','nonnegative','2d','size',[Ni,u]});
        end
        if isdef.fitSPI, opt.fitSPI = false; end
        isdef.SPI = false;
    end
    
    if opt.normalizeSPI, opt.SPI = opt.SPI./sum(opt.SPI,2).*opt.nSPI'; end
    
    if any(abs(sum(opt.SPI,2)./opt.nSPI'-1) > 0.2)
        warning('SPI seems inconsistent with nSPI');
    end

    if parsingpending
        parsestruct(data,FLD.mat{1},'opt',FLD.mat{2},'real','2d','size',[NaN,Ni]);
    end
    
    if opt.TmWeight > 0
        
        % Assume temperature sensor, if available, is at the center of the table
        if isempty(opt.TmRow), opt.TmRow = ceil((opt.rows+1)/2); end
        
        if isscalar(opt.TmRow), opt.TmRow = full(interpmatrix(opt.TmRow,(1:opt.rows)')); end
            
        if u > opt.rows && numel(opt.TmRow) ~= u
            try
                jj = find(opt.TmRow);
                jj = arrayfun(@(k) find(all(circshift(eye(opt.rows,1),k-1)' == opt.strtypes,2)),jj);
                opt.TmRow = full(sparse(1,jj,nonzeros(opt.TmRow),1,size(opt.strtypes,1)));
            catch
                error('Failed to project TmWeight in terms of string types');
            end
        end        
        validateattributes(opt.TmRow,'numeric',{'real','2d','finite','nonnegative','size',[1,u]});
    end
    
    % 1 block-per-row for vertical modules, 3 blocks per row for horizontal modules
    if isdef.BPR, opt.BPR = 1 + 2*(opt.mdims(1) > opt.mdims(2)); end
    
    if opt.dryrun && nargout <= 1, varargout = {opt}; return; end

    L = opt.rows*opt.mdims(2) + (opt.rows-1)*opt.mgap;  % table width (along module surface)
    assert(opt.pitch > L.*cosd(opt.tilt),'Expecting pitch > %0.2f',L.*cosd(opt.tilt));
    
    A = prod(opt.mdims); % module area
    assert(opt.PmpRef/(A*1000) < 0.5,'Check PmpRef and mdims, efficiency > 50%');
    
    MPI = opt.MPS*opt.nSPI; % modules per inverter
    
    if opt.verbose
       SFA = infinite_rows(opt.pitch,L,opt.tilt,1);
       fprintf('Table width: %0.2f m\nShade-free angle: %0.1f°\n',L,SFA);
    end
    
    if opt.fitinverter

        if max(data.Pac(:)) > opt.PacMax
            warning('Monitoring values with up to %0.1f Pac0',max(data.Pac(:))/opt.PacMax);
        end

        f = data.Pac > 0.1*opt.PacMax & data.Pac < 0.95*opt.PacMax; % avoid clipping
        if nnz(f) > 1e5, f(f) = rand(nnz(f),1) < 1e5/nnz(f); end
        mdl = fit(data.Pdc(f),data.Pac(f),'(x-Ps0)*EffMax','robust','bisquare','StartPoint',[0 0.95]);
        
        if opt.verbose
            fprintf('\nInverter efficiency: %0.2f%%\n',mdl.EffMax*100);
            fprintf('Inverter offset: %0.1f W\n',mdl.Ps0);
        end

        opt.Inverter(1).EffMax = mdl.EffMax;
        opt.Inverter.Ps0 = mdl.Ps0;
        opt.Inverter.Pdc0 = opt.PacMax/mdl.EffMax;
        % opt.Inverter.fPV  = @(Pdc,Vdc) max(0,Pdc - mdl.Ps0)*mdl.EffMax;

        f = data.Pac > 0.001;
        if ~isempty(data.Vdc)
            opt.Inverter.MPPTHi = max(data.Vdc(f));
            opt.Inverter.MPPTLow = min(data.Vdc(f));
        else
            opt.Inverter.MPPTHi = Inf;
            opt.Inverter.MPPTLow = 0;
        end
        if ~isempty(data.Idc)
            opt.Inverter.IdcMax = max(data.Idc(f));
        else
            opt.Inverter.IdcMax = Inf;
        end

        % INV.Pnt
        % Inv.Paux
        % INV.TPLim = @(Ta) temperature clippint
    end
    
    if ~isempty(data.tpw) && ~isempty(data.AMa)
        FS = pvl_FSspeccorr(data.tpw/10, data.AMa,opt.pvModType);
    else
        warning('Missing tpw and AMa to calculate spectral correction');
        FS = 1;
    end
    
    fIAM = checkIAM('martinruiz','maxloss',opt.IAM_maxloss);
    
    % Infinite line (2D) isotropic  shading
    [Ge,ShF,DTI,BTI0,GTIm] = infinite_rows(opt.pitch,L,opt.tilt,opt.rows,data,opt.azimuth,...
        'GTI',opt.knownGTI,'albedo',opt.albedo,'fIAM',fIAM,'CSR',opt.CSR);

    M = double(opt.strtypes./sum(opt.strtypes,2));
    
    G0 = FS.*(DTI + BTI0)*M';
    Ge = FS.*Ge*M';
    DTI = FS.*DTI*M';
    
    % Approx. Electrical shading losses
    Nsb = ceil(ShF*opt.strtypes'*opt.BPR);  % shaded blocks per string
    Ntb = opt.BPR*sum(opt.strtypes,2)';     % total blocks per string
    
    if strcmp(opt.mismatch,'deline')
        Fe = inverter_mismatch(Nsb,Ntb,DTI./G0,opt); % [t,u,Ni] (one plane for each inverter)
    else
        Fe = mismatch(Nsb,Ntb,max(0,DTI./G0),'model',opt.mismatch);
    end

    if ~isempty(data.Pdc)
        data.Pdc(data.Pdc <= 0) = NaN;
        data.Pdc = data.Pdc./(MPI.*opt.PmpRef); % W/Wp
        
        PR = data.Pdc./GTIm*1000;
        data.eta = PR*(opt.PmpRef/1000)/A;
    else
        data.eta = [];
    end

    if ~isempty(data.Pac)
        data.Pac(data.Pac <= 0) = NaN;
        data.Pac = data.Pac./opt.PacMax;
    end

    ok = data.sunel > opt.min_sunel & GTIm > opt.min_GTI; 
    if ~isempty(data.Pdc)
        ok = ok & data.eta < opt.max_eta & ...
                  data.Pdc > opt.min_Pdc;
    end
    if ~isempty(data.Pac)
        ok = ok & data.Pac > opt.min_Pac & ...
                  data.Pac < opt.max_Pac;
    end
        
    switch mode
    case 'fit'
                
        % Divide data into train-test sets, filter out samples with clipping
        [~,~,d] = unique(floor(solarposition.doy(data.t)-0.5));
        ff = ok & mod(d,3) == 0 & data.Pdc < opt.clip_Pdc.*opt.Inverter.Pdc0./(MPI.*opt.PmpRef);
        if isempty(data.Vdc)
            if isfield(opt.Inverter,'MPPTHi'), ff = ff & data.Vdc <= opt.Inverter.MPPTHi; end
            if isfield(opt.Inverter,'MPPTLow'), ff = ff & data.Vdc >= opt.Inverter.MPPTLow; end
        end
        
        bestRMSE = Inf;
        bestSPI = opt.SPI;
        
        for iter = 1:10
            f = all(ff,2);

            if iter == 1 && opt.fitSPI
            % For the first run, stick to (nearly) non-shaded samples
                f = f & all(Fe > 0.95,2:3);
            end

            % Fit model parameters
            [p,resnorm] = model_fit(Ge(f,:),data.Ta(f),data.vw(f),Fe(f,:,:),data.Pdc(f,:),data.Tm(f,:),opt);
            
            % % DEBUG
            % [eta_mdl,Tm_mdl,Pdc_mdl] = model_run(p,Ge(f,:),data.Ta(f),data.vw(f),Fe(f,:),opt);
            % GUIfigure('PVNOR_fit'); clf(); prettyscatter(Pdc_mdl,data.Pdc(f,:),[],'Pdc [W/Wp]',@mean,'%');                
            % keyboard();

            if resnorm > bestRMSE
                opt.SPI = bestSPI;
                break;
            else
                if iter > 1, bestRMSE = resnorm; end
                bestSPI = opt.SPI;
            end
            
            % Fit inverter row mix matrix
            if ~opt.fitSPI, break; end
            opt.SPI = fitSPI(p,Ge(f,:),data.Ta(f),data.vw(f),Fe(f,:,:),data.Pdc(f,:),ff(f,:),opt);

            if size(Fe,3) > 1
                Fe = inverter_mismatch(Nsb,Ntb,DTI./G0,opt);
            end
        end
        
        opt.params = p;
        
        ff = mod(d,3) ~= 0 & ok;
    case 'run'
        ff = ok;
        p = opt.params;
        
    otherwise, error('this should not happen');
    end

    f = any(ff,2);
    [~,Tm_mdl,Pdc_mdl,Pac_mdl] = model_run(p,Ge(f,:),data.Ta(f),data.vw(f),Fe(f,:,:),opt);

    fprintf('\n%s\n',modelname(p,opt.model,1));
    
    [MBE,RMS,MAD] = printstats({Tm_mdl,Pdc_mdl,Pac_mdl},{data.Tm(f,:),data.Pdc(f,:),data.Pac(f,:)},...
                               {'Tm','Pdc','Pac'});

    if opt.plot
        n = round(sqrt(9*Ni/16));
        m = ceil(Ni/n);
        GUIfigure('PVsimple','PV_simple',sprintf('%d:%d',m,n)); clf();

        if opt.normalizeSPI
            tags = arrayfun(@(j,l) sprintf('INV %d',j,l),1:Ni,'unif',0);
        else
            loss = (1-sum(opt.SPI,2)./opt.nSPI')'*100;
            tags = arrayfun(@(j,l) sprintf('INV %d (%0.1f%% loss)',j,l),1:Ni,loss,'unif',0);
        end
        multiscatter(Pac_mdl,data.Pac(f,:),n,m,'titles',tags,'densityplot',true,...
            'xlabel','P_{ac,mdl} [kW/kW_0]','ylabel','P_{ac,mes} [kW/kW_0]','oversize',[1.1,1.1]);
    end

    f = revertfilter(f,useful);
    varargout = {opt,MBE,RMS,MAD,f};
end
    
function [eta_mdl,Tm_mdl] = mattei(p,Ge,Ta,vw)

    c = num2cell(p); 
    [uconst,uwind,alpha,eta_ref,beta,gamma] = deal(c{:});
    
    T_STC = 25;
    
    Upv = uconst + uwind*vw;
    
    eta_mdl = eta_ref.*(1 + gamma*log10(Ge/1000) + beta*T_STC);
    
    Tm_mdl = (Upv.*Ta + Ge.*(alpha - eta_mdl))./(Upv - beta*eta_ref*Ge);
    eta_mdl = eta_mdl - beta*eta_ref.*Tm_mdl;
    
    % eta_mdl(eta_mdl > 0.5) = NaN;
end

function [eta_mdl,Tm_mdl] = faiman(p,Ge,Ta,vw,eff)

    k = p(4); % mutiplying factor for soiling, mismatch, ohmic, etc.
    p = cell2struct(num2cell(p(1:3))',{'Uconst','Uwind','absort'});

    Tm_mdl = celltempUvalues(Ge,Ta,vw,p,@(g,t) k*eff(g,t),1e-3);
    eta_mdl = k*eff(Ge,Tm_mdl);
end

function [eta_mdl,Tm_mdl,Pdc_mdl,Pac_mdl] = model_run(p,Ge,Ta,vw,fe,opt)

    switch opt.model
    case 'mattei'
        % validateattributes(p0,'numeric',{'vector','real','numel',6});
        [eta_mdl,Tm_mdl] = mattei(p,Ge,Ta,vw);
        Vmp = [];
    case 'odm+faiman'
        % validateattributes(p0,'numeric',{'vector','real','numel',4});
        [eta_mdl,Tm_mdl] = faiman(p,Ge,Ta,vw,@(g,t) opt.Module.MPPefficiency(g,t));
        
        if nargout > 3
            [~,Vmp] = opt.Module.getMPP(Ge,Tm_mdl,'interp');
            Vmp = mean(Vmp*opt.MPS,2);
        end
        
    otherwise, error('this should not happen');
    end

    eta_mdl = eta_mdl.*fe;

    if nargout > 2
    % Pdc = (eta·Ge)·SPI'·MPS·A / (nSPI·MPS·PmpRef) = (eta·Ge)·SPI'·A/(nSPI·PmpRef) [W/Wp]
        K = prod(opt.mdims)./(opt.nSPI*opt.PmpRef);
        if size(eta_mdl,3) == 1
            Pdc_mdl = (eta_mdl.*Ge)*opt.SPI'.*K;
        else
            Pdc_mdl = sum(permute(eta_mdl.*Ge,[1,3,2]).*permute(opt.SPI,[3,1,2]),3).*K;
        end      
    end
    if nargout > 3
        WpPI = opt.PmpRef*opt.MPS*opt.nSPI;
        Pac_mdl = ONDeval(opt.Inverter,Pdc_mdl.*WpPI,Vmp,Ta)/opt.PacMax;
    end
end

function p0 = modeldefaults(model)
    switch model
    case 'mattei'
        [~,~,p0] = celltemp_mattei(1000,25,1);
        p0 = struct2cell(p0);
        p0 = cat(2,p0{:});
        % p0 = [24.1,2.9,0.81,0.2,0.005,0]; 
    case 'odm+faiman'
        p0 = [24.1,2.9,0.81,1];
        
    otherwise, error('this should not happen');
    end
end

function msg = modelname(p,model,long)

    if nargin < 3, long = false; end
    
    p0 = modeldefaults(model);
    
    if isequal(p,p0), msg = 'default'; else, msg = 'custom'; end
    msg = [upper(model(1)),model(2:end) ' - ' msg];
    
    if ~long, return; end
    
    switch model
    case 'mattei'
        msg = sprintf('%s (Upv = %0.1f + %0.1f v, a = %0.2f, eta0 = %0.2f, b = %0.3f, g = %0.2f)',msg,p);
    case 'odm+faiman'
        msg = sprintf('%s (Upv = %0.1f + %0.1f v, a = %0.2f, k = %0.2f',msg,p);
        
    otherwise, error('this should not happen');
    end
end

function Fe = inverter_mismatch(Nsb,Ntb,rd,opt)
% for mismatch = 'deline', estimates fraction of shaded strings X based on the current SPI
    
    Ni = size(opt.SPI,1);
    X = zeros([size(Nsb),Ni]);
    for j = 1:Ni
        X(:,:,j) = sum((Nsb >= permute(Nsb,[1,3,2])).*shiftdim(opt.SPI(j,:),-1),3)./opt.nSPI(j);
    end
    assert(all(X < 1 + eps(Ni) & X > -eps(Ni),'all'));
    X = max(0,min(1,X));

    Fe = mismatch(Nsb,Ntb,rd,'model',opt.mismatch,'X',X,opt.mismatch_opt{:});
end

function [p,resnorm] = model_fit(Ge,Ta,vw,fe,Pdc_mes,Tm_mes,opt)

    % if nargin < 7, Tm_mes = []; end
    % compatiblesize(Ge,Ta,vw,fe,Pdc_mes);

    [ni,u] = size(opt.SPI);
    assert(size(Ge,2) == u);
    assert(size(fe,2) == u);
    if size(fe,3) > 1, assert(size(fe,3) == ni); end
    
    if opt.verbose
        fprintf('\nFitting %s (%d samples) ...\n',opt.model,size(Ge,1));
    end

    Pdc_mdl = [];
    Tm_mdl = [];

    p0 = opt.params;
    
    if contains(opt.fitmethod,'g')
        n = round(sqrt(size(Ge,1)));
        [~,G] = sort(Ge);
        G = ceil(G/n);
    end

    if opt.verbose
        model_err(p0); % sets eta_mdl, Tm_mdl, Pdc_mdl
        
        fprintf('\n%s\n',modelname(p0,opt.model,1));
        printstats({Pdc_mdl,Tm_mdl},{Pdc_mes,Tm_mes},{'Pdc','Tm'});
    end

    optopt = optimoptions(@lsqnonlin,'MaxFunEvals',2e3,'Display','none',...
        'FiniteDifferenceStepSize',1e-4,'DiffMinChange',1e-6);
    
    switch opt.model
    case 'mattei'
        % [uconst,uwind,alpha,eta_ref,beta,gamma]
        lb = [0 0 0 0 0 -1];    % lower bound
        ub = [30,10,1,1,0.1 1]; % upper bound
        
        optopt.TypicalX = [20,5,1,0.1,0.01 1];
    case 'odm+faiman'
        % [uconst,uwind,alpha,absort]
        lb = [0 0 0 0];
        ub = [30,10,1,1];
                
        optopt.TypicalX = [20,5,1,1];
        optopt.Display = 'iter';
        % optopt.DiffMinChange = 1e-2;
        % optopt.StepTolerance = 1e-2;
    otherwise, error('this should not happen');
    end

    [p,resnorm,~,exitflag]  = lsqnonlin(@model_err,p0,lb,ub,optopt);
    assert(exitflag > 0);

    % % Run lsqnonlin(fun,p0,lb,ub), starting from 10 different seeds
    % problem = createOptimProblem('lsqnonlin','x0',p0,'lb',lb,'ub',ub,'objective',@mattei);
    % ms = MultiStart();
    % p = run(ms,problem,10)

    if opt.verbose
        model_err(p);
        
        fprintf('\n%s\n',modelname(p,opt.model,1));
        printstats({Pdc_mdl,Tm_mdl},{Pdc_mes,Tm_mes},{'Pdc','Tm'});
    end

    function err = model_err(p)

        [~,Tm_mdl,Pdc_mdl] = model_run(p,Ge,Ta,vw,fe,opt);
        
        % Pdc_mdl = getPdc(eta_mdl,Ge.*fe);
        err = double(Pdc_mdl - Pdc_mes);
        % err = double((eta_mdl - eta_mes).*Ge);
        % err = double(eta_mdl.*fe - eta_mes);
        
        if ~isempty(Tm_mes) && opt.TmWeight > 0
            Tm_mdl = Tm_mdl*opt.TmRow';
            err(:,end+1) = double((Tm_mes - Tm_mdl)*opt.TmWeight.*(Ge*opt.TmRow')/1000);
        end
        
        switch opt.fitmethod
        case 'lsq'
        case 'lar', err = sign(err).*sqrt(abs(err));
        case 'lgsq', err = splitapply(@mean,err,G);
                
        otherwise, error('this should not happen');
        end
    end
end

function SPI = fitSPI(p,Ge,Ta,vw,Fe,Pdc_mes,ok,opt)

    % Pdc = (eta·Ge)·SPI·MPS·A / (nSPI·MPS·PmpRef) = eta·Ge·(SPI/nSPI)·(A/PmpRef) = eta·Ge·w·nA 
    nA = prod(opt.mdims)/opt.PmpRef; 
            
    Ni = numel(opt.nSPI);
    u = size(opt.strtypes,1);
    SPI = zeros(Ni,u);
    
    if opt.verbose
        nok = sum(ok,1);
        fprintf('\nFitting SPI (%d-%d samples per inverter) ...\n',min(nok),max(nok));
    end
    
    m = size(Fe,2);
    assert(m == u && size(Ge,2) == m && size(Pdc_mes,2) == Ni);
    if size(Fe,3) > 1
        assert(size(Fe,3) == Ni);
    else
        Fe = repmat(Fe,1,1,Ni);
    end
    
    for j = 1:Ni
        f = ok(:,j);
        eta_mdl = model_run(p,Ge(f,:),Ta(f),vw(f),1,opt);
        
        Pdc_row = eta_mdl.*Ge(f,:).*Fe(f,:,j)*nA; % normalized Pdc @ nominal nSPI
        w = fitrowmix(Pdc_row,Pdc_mes(:,j));
        SPI(j,:) = opt.nSPI(j)*w';
        % SPI = SPI./sum(SPI,2);
    end
    
    % Pdc_mdl = (eta_mdl.*Ge)*SPI'./(nSPI)*nA

    function w = fitrowmix(X,y)

        [w_nn,~,~,exitflag] = lsqnonneg(double(X),y);
        if exitflag <= 0, w_nn = ones(1,m)/m; end

        if contains(opt.fitmethod,'g')
            n = round(sqrt(numel(y)));
            [~,G] = sort(y);
            G = ceil(G/n);
        end

        optopt = optimoptions(@lsqnonlin,'MaxFunEvals',1e5,'Display','none',...
            'FiniteDifferenceStepSize',1e-4,'StepTolerance',1e-3);
        [w,~,~,exitflag] = lsqnonlin(@errfn,w_nn,zeros(m,1),[],optopt);
        assert(exitflag > 0);

        if opt.normalizeSPI || sum(w) > 1, w = w/sum(w); end

        function err = errfn(w)
            
            s = sum(w);
            if opt.normalizeSPI || s > 1, w = w/sum(w); end
            err = double(X*w - y);

            switch opt.fitmethod
            case 'lsq'
            case 'lar', err = sign(err).*sqrt(abs(err));
            case 'lgsq', err = splitapply(@mean,err,G);
                    
            otherwise, error('this should not happen');
            end
        end
    end
end

function varargout = printstats(mdl,mes,typ)
    
    TYPE = {'eta','Tm','Pdc','Pac'};
    NORM = {1 1 1e-3 1e-3};
    UNITS = {'%','K','W/kWp','W/kW0'};
    
    [lbl,idx] = parselist(typ,TYPE);
    norm = NORM(idx);
    units = UNITS(idx);
    
    for j = 1:numel(mdl)
        if isempty(mes{j}), continue; end
        [mdl{j},mes{j}] = compatiblesize(mdl{j},mes{j});
        [msg,mbe,rms,mad] = errorstats(mdl{j}(:),mes{j}(:),norm{j},'units',units{j});
        fprintf('\t%s: %s\n',lbl{j},msg);
        
        if nargout > 0
            MBE.(lbl{j}) = mbe;
            RMS.(lbl{j}) = rms;
            MAD.(lbl{j}) = mad;
        end
    end
    if nargout > 0
        varargout = {MBE,RMS,MAD};
    end
end
