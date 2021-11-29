function Inverter = ONDread(filename)
% INV = ONDREAD(FILE) - Reads an *.OND PVsyst-inverter model file, or a text file (*.ondcp) 
%   containing an exported model (semicolon delimited values, as in PVsyst 6.7.5).
%   The details of the format (order of fields, and aliases for *.OND files) can be found in
%   config file ONDFORMAT.TXT.
%
%   PROVISIONAL: Only works with voltage-dependent inverters(i.e. explicit curves).
%
%   Generates a structure INV to be used with ONDVAL.

    if nargin < 1, filename = '*.ondcp'; end
    filename = pickfile(filename,'Pick an ODM file');
    
    DEF.Pnt = 0;
    DEF.Paux = 0;
    
    if str2double(regexprep(version(),'(\d+\.\d+).*','$1')) <= 9.2
        METHOD_EFF = 'pchip';
    else
        METHOD_EFF = 'makima';
    end
    
    if ~isempty(regexpi(filename,'.*\.OND$'))
        Inverter = ONDread_legacy(filename);
    else
        Inverter = ONDread_copied(filename);
    end
    % assert(~isempty(Inverter.name),'Unrecognized inverter file format');
    
    if Inverter.nMPPT > 1 
    % PROVISIONAL: split inverter into nMPPT independent 'virtual' inverters
    
        warning('Splitting inverter into %d independent MPPTs',Inverter.nMPPT)

        % 'info.MasterSlave','info.NbMSInterne','info.NbInputs','nMPPT'
    
        MPPT_SCALABLE = {'info.INomAC','IacMax','info.INomDC','IdcMax',...
                 'info.PNomAC','PacMax','info.PNomDC','Pdc0','Ps0','Paux','PauxThresh','Pnt'};
             
        for f = MPPT_SCALABLE
            if ~isnestedfield(Inverter,f{1}), continue;
            else
                v = getnestedfield(Inverter,f{1});
                Inverter = setnestedfield(Inverter,f{1},v/Inverter.nMPPT);
            end
        end
        if isfield(Inverter,'TPLim') && size(Inverter.TPLim,2) > 1
            Inverter.TPLim(:,2) = Inverter.TPLim(:,2)/Inverter.nMPPT;
        end
        if isfield(Inverter,'refcurve')
            Inverter.refcurve(:,1) = Inverter.refcurve(:,1)/Inverter.nMPPT;
        end
        for k = 1:numel(Inverter.curves)
            Inverter.curves{k}(:,1) = Inverter.curves{k}(:,1)/Inverter.nMPPT;
        end
        Inverter.info.Model = sprintf('1/%d %s',Inverter.nMPPT,Inverter.info.Model);
        Inverter.info.SCALED = sprintf('1/%d MPPT FRACTION OF ORIGINAL INVERTER!',Inverter.nMPPT);
        Inverter.nMPPT = 1;
        %'info.FileName','info.DataSource',
    end
    
    Inverter.name = Inverter.info.Model;
    assert(Inverter.vdependent,'Sorry, only working with 3-voltage curves!');

    % Temperature-clipping interpolant between points:
    %       TPLim = [     0,  TPMax,  TPNom, TPLim1, TPLimAbs ;
    %                PacMax, PacMax, PNomAC,  PLim1, PLimAbs ]'
    if Inverter.hasTPLim
       TPLim = Inverter.TPLim;
       TPLim(1,1:2) = [0,Inverter.PacMax];
       TPLim(2:3,2) = [Inverter.PacMax;Inverter.info.PNomAC];
       TPLim(6:7,1) = TPLim(5)+[1e-6,10];
       Inverter.TPLim = griddedInterpolant(TPLim(:,1),TPLim(:,2),'linear','nearest');
    else
       Inverter.TPLim = @(x) ones(size(x))*Inverter.PacMax;
    end

    if ~isfield(Inverter,'IdcMax') && all(isfield(Inverter,{'VdcMin','EffMaxV'}))
        Inverter.IdcMax = min(Inverter.IdcMax,Inverter.PacMax/min(Inverter.EffMaxV)/Inverter.VdcMin);
    end
    
    Inverter = completestruct(Inverter,DEF);
    
    % Convert Inverter-Curves to Input-Power & Voltage interpolant
    
    % Resample curves for smooth gridded-interpolant
    pts = Inverter.curves;
    pq = [0;linspace(Inverter.Ps0,Inverter.Pdc0,1000)']; % note starting zero
    zq = zeros(numel(pq),numel(pts));
    for j = 1:numel(pts)
        
        % Ensure the curve starts at zero-power and can reach PacMax
        if pts{j}(1,2) > 0 && pts{j}(1,1) > Inverter.Ps0
            pts{j} = cat(1,[Inverter.Ps0,0],pts{j});
        end
        if pts{j}(end,1) < Inverter.PacMax
            x0 = interp1(pts{j}(:,2).*pts{j}(:,1),pts{j}(:,1),Inverter.PacMax,'linear','extrap');
            pts{j} = cat(1,pts{j},[x0,Inverter.PacMax/x0]);
        end
        
        % Use a weighted mixture of Akima- (on efficiency) and Linear- (on Power) Interpolation
        % setting w ~ eff ensures no overshots at high-efficiency, and a smooth low-power curve
        E = griddedInterpolant(pts{j}(:,1),pts{j}(:,2),METHOD_EFF,'nearest');
        P = griddedInterpolant(pts{j}(:,1),pts{j}(:,2).*pts{j}(:,1),'linear','nearest');
        w = P(pq); 
        w = w./max(w);
        zq(:,j) = w.*(E(pq).*pq) + (1-w).*P(pq);
    end
    zq(1,:) = 0;
    
    if size(zq,2) > 1
        [Vc,idx] = sort(Inverter.Vc(:));
        zq = zq(:,idx);
        
        % Add nearest-neighbor curves at MPPT-boundaries (avoid need for extrapolation)
        if Inverter.MPPTLow < Vc(1)
            Vc = [Inverter.MPPTLow; Vc];
            zq = [zq(:,1),zq];
        end
        if Inverter.MPPTHi > Vc(end)
            Vc = [Vc; Inverter.MPPTHi];
            zq = [zq,zq(:,end)];
        end
        [pq,vq] = ndgrid(pq,Vc);
        
        Inverter.fPV = griddedInterpolant(pq,vq,zq,'linear');
    else
        fPV = griddedInterpolant(pq,zq,'linear');
        Inverter.fPV = @(p,v) fPV(p);
    end

    Inverter.fOut = @(p,v,t) ONDeval(Inverter,p,v,t);
end

function Inverter = ONDread_copied(filename)
% Read a semicolon-delimited list of propperties, as exported to the clipboard by PVSYST (6.7.5)
% The meaning (i.e. order) of the propperties is defined in config. file "ONDformat.txt"

    % ONDformat.txt should contain at least three columns: varname, type, scale,...  
    fileID = fopen('ONDformat.txt');
    if fileID < 0, error('ONDread: Could not find ONDformat.txt'); end
    filecloser = onCleanup(@() fclose(fileID));
    fmt = textscan(fileID,'%s %s %f %s','delimiter',';','multipledelimsasone',false,'commentstyle','#');
    delete(filecloser); % fclose(fileID);
    fmt(4) = [];
    fmt{3}(isnan(fmt{3})) = 1;
    
    % Blocks of voltage-dependent curves, starting with Vc(n); ... are optional
    oksizes = [find(cellfun(@(s) ~isempty(regexp(s,'^Vc\(\d\).*','once')),fmt{1}))-1;numel(fmt{1})];
    
	fileID = fopen(filename);
    if fileID < 0, error('ONDread: Could not open file'); end
    filecloser = onCleanup(@() fclose(fileID));
    data = textscan(fileID,'%s','delimiter',';','multipledelimsasone',false);
    delete(filecloser); % fclose(fileID);
    data = data{1};
    
    assert(any(numel(data) == oksizes),['File does not match known format (PVSyst 6.7.5), ',...
        'Review/modify config. file "ONDformat.txt" to reflect OND spec. changes']);
    
    floats = strcmp(fmt{2}(1:numel(data)),'float');
    decimalcomma = nnz(contains(data(floats),',')) > nnz(contains(data(floats),'.'));

    Inverter = struct('name','');
    for j = 1:numel(data)
        switch fmt{2}{j}
            case 'bool', val = ~isempty(data{j});
            case 'yn'
                if isempty(data{j}), val = NaN;
                else, val = lower(data{j}(1))=='y';
                end
            case {'enum','str'}, val = data{j};
            case 'float'
                if decimalcomma, data{j} = strrep(data{j},',','.'); end
                val = str2double(data{j})*fmt{3}(j);
            case 'int', val = str2double(data{j});
        end
        if ~isnan(val)
            eval(['Inverter.' fmt{1}{j} '= val;']);
        end
    end
    
    % Assign name if things went smooth
    if isnestedfield(Inverter,'info.FileName') && ~isempty(Inverter.info.FileName)
        Inverter.name = strrep(Inverter.info.FileName,'.OND','');
    elseif isnestedfield(Inverter,'info.Model') && ~isempty(Inverter.info.Model)
        Inverter.name = Inverter.info.Model;
    end
end

function Inverter = ONDread_legacy(filename)
% Reads an *.OND PVsyst-inverter model file (KML-feel), e.g:
%
%     PVObject_=pvGInverter
%     ...
%     PVObject_Commercial=pvCommercial
%     ...
%     End of PVObject pvCommercial
%     Transfo=Without
%     ...
%     Converter=TConverter
%         PNomConv=25.000
%         PMaxOUT=25.000
%         ...
% 
%         ProfilPIOV1=TCubicProfile
%           ...
%         End of TCubicProfile
%     End of TConverter
%     NbInputs=8
%     ...
%     End of PVObject pvGInverter
%
% Known field names are contained in the 4th column of "ONDformat.txt" (ALIAS)

    if nargin < 1, filename = '*.OND'; end
    filename = pickfile(filename,'Pick an ODM file');
    
	fileID = fopen(filename);
	if fileID < 0, error('ONDread: Could not open file'); end
    
    % Read info header lines
    data = textscan(fileID,'%[^\n\r]');
    fclose(fileID);
    
    Inverter = struct('name','');
    
    % Read curves, wrapped in {ProfilPIO=TCubicProfile,... End of TCubicProfile} blocks
    [Inverter.curves,data] = readprofiles(data{1});
    [idx,vals] = getindex(data);
    
    approx = @(s) regexprep(lower(strtrim(s)),'[\W_]','');
    
    % Read [varname, type, scale, alias] from ONDformat.txt
    fileID = fopen('ONDformat.txt');
    if fileID < 0, error('ONDread: Could not find ONDformat.txt'); end
    filecloser = onCleanup(@() fclose(fileID));
    fmt = textscan(fileID,'%s %s %f %s','delimiter',';','multipledelimsasone',false,'commentstyle','#');
    delete(filecloser); % fclose(fileID);
    fmt{3}(isnan(fmt{3})) = 1;
    
    % Generate an extended list of aliases
    alias = cellfun(@(s) split(s,'/'),fmt{4},'unif',0);
    ia = arrayfun(@(j) repmat(j,numel(alias{j}),1),1:numel(alias),'unif',0);
    alias = cat(1,alias{:});
    ia = cat(1,ia{:});
    
    % Find matches for variable aliases
    [known,ic] = ismember(approx(idx),approx(alias));
    ic(known) = ia(ic(known));
    
    if any(known)
        floats = known;
        floats(known) = strcmp(fmt{2}(ic(known)),'float');
        decimalcomma = nnz(contains(vals(floats),',')) > nnz(contains(vals(floats),'.'));

        % Create structure from recognized field names
        for j = find(known)'
            switch fmt{2}{ic(j)}
                case 'bool', val = ~isempty(vals{j});
                case 'yn'
                    if isempty(vals{j}), val = NaN;
                    else, val = lower(vals{j}(1))=='y';
                    end
                case {'enum','str'}, val = vals{j};
                case 'float'
                    if decimalcomma, vals{j} = strrep(vals{j},',','.'); end
                    val = str2double(vals{j})*fmt{3}(ic(j));
                case 'int', val = str2double(vals{j});
            end
            if ~isnan(val)
                eval(['Inverter.' fmt{1}{ic(j)} '= val;']);
            end
        end
    end

    Inverter.Vc = findfield(idx,vals,'VNomEff',[],'%f','Delimiter',',','CollectOutput',1);
    Inverter.EffMaxV = findfield(idx,vals,'EfficMaxV',[],'%f','Delimiter',',','CollectOutput',1)'/100;
    Inverter.EffEuroV = findfield(idx,vals,'EfficEuroV',[],'%f','Delimiter',',','CollectOutput',1)'/100;
    
    Inverter.hasTPLim = isfield(Inverter,'TPLim');
    Inverter.vdependent = numel(Inverter.Vc) > 1;
end

function [idx,vals] = getindex(data)
% PROVISIONAL: should identify nested features
%              should handle things like 'Remarks, count = 4', 'End of...', etc.

    Nf = numel(data);
    idx = cell(Nf,1);
    vals = cell(Nf,1);
    empty = false(Nf,1);
    
    for j = 1:Nf
        k = strfind(data{j},'=');
        empty(j) = isempty(k);
        if ~empty(j)
           idx{j} = data{j}(1:k-1); 
           vals{j} = data{j}(k+1:end); 
        end
    end
    
    idx = idx(~empty);
    vals = vals(~empty);
end

function [pts,rest] = readprofiles(data)
% Search OND file for patterns:
%
%     ProfilPIO=TCubicProfile
%       NPtsMax=11
%       NPtsEff=8
%       LastCompile=$0011
%       Mode=1
%       Point_1=19800.0,0.0
%       ...
%       Point_8=1980000.0,1955871.2
%       ...
%       Point_11=0.0,0.0
%     End of TCubicProfile
%
% ... and return them as a cell-array of [N,2] curves (PDC vs Eff), along with the ~used data.

    used = false(numel(data),1);
    
    ends = find(strcmp(data,'End of TCubicProfile'));
    starts = find(~cellfun(@isempty,regexp(data,'.*ProfilPIO+V?(\d)?=TCubicProfile+')));

    pts = cell(numel(starts),1);
    for j = 1:numel(starts)
        used(starts(j):ends(j)) = true;
        [idx,vals] = getindex(data(starts(j)+1:ends(j)-1));
        Np = findfield(idx,vals,'NptsEff');
        if findfield(idx,vals,'Mode')~= 1, warning('Curves might not make sense'); end

        pts{j} = zeros(Np,2);
        for i = 1:Np
            pts{j}(i,:) = findfield(idx,vals,sprintf('Point_%d',i),[],'%f,%f','CollectOutput',1);
        end
        pts{j}(:,2) = pts{j}(:,2)./pts{j}(:,1);
    end
    
    rest = data(~used);
end

function [val,k] = findfield(idx,vals,fname,default,formatstring,varargin)
% Find the field 'fname' in idx, and return its corresponding value (and optionally line number k)
% If the field is not found, default will be returned, if it is not available, an error will be raised
% If formatstring is provided, it will be passed to textscan(str,formatstring), along with varargin

    k = find(strcmpi(fname,idx));
    if isempty(k)
       if nargin > 3, val = default; k = 0; return;
       else
           error('ONDread:findfield','field ''%s'' not found, no default provided',fname); 
       end
    end

    if nargin < 5, formatstring = '%f'; end

    val = textscan(vals{k},formatstring,varargin{:});
    if numel(val) == 1, val = val{1}; end

end