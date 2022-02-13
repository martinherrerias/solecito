function adjust2specs(prjname,type)
% ADJUST2SPECS(PRJNAME,TYPE) - Generate modified copies of files mount.TXT, positions.TXT, 
%   module.TXT, and Array Definition.TXT generated by the preprocessor, that (more closely) 
%   follow the specifications defined in May 2016.
%
%   In other words, it does the job of other people, too slow to follow simple instructions.
%
% Updated 9.8.2019, tested with 0a output.
% TODO: check/adjust for 1aC preprocessor output.

    nicestring = @(s) matlab.lang.makeValidName(regexprep(regexprep(strrep(s,' ','_'),'_+','_'),'(.*?)_?$','$1'));

    layoutfile = pickfile('*positions*.TXT');
    PATH = fileparts(layoutfile); cd(PATH); 
    modfile = pickfile('*module*.TXT');
    mountfile = pickfile('*mount*.TXT');
    arrdeffile = pickfile('*Definition*.TXT');
    
    LEGEND = 'Adjusted automatically to specification: check for compliance/completeness!';
    PREFIX = nicestring(prjname);
    
%% Fix physical module *.mpoly file:

    % #file_format: PREPROCESSOR_V2.0
    % #module_id: [nonsense (discarded)]
    % #Module TSM-PD14 72 Cells [taken as base for true module_id]
    % #
    % #outline_polygone: 0;0;0.976;0;0.976;1.948;0;1.948 [sic]
    % #divisions(n;string_id;x1;y1;x2;y2;..)
    % 1;1;0.005;0.005;0.161;0.005;0.161;0.161;0.005;0.161 [zero-based indices]
    
    module = readtxtfile(modfile,{'Module','module_id';'outline_polygone:','outline_polygon'},'colon',' ','CollectOutput',0);
    module.comments = {LEGEND;module.params.module_id;moduleinfo(module)};
    module.params.module_id = nicestring(module.params.module_id);
    module.headers = {'divisions(n;string_id;x1;y1;x2;y2;...)'};
    module.data{1} = module.data{1} + 1; % fix zero-based indices
    module.data{2} = module.data{2} + 1;
    modfile = fullfile(PATH,[module.params.module_id '.mpoly']);
    % format = ['%02d;%01d;' strjoin(repmat({'%0.3f'},1,numel(module.data)-2),';')];
    writetxtfile(modfile,module);
    
%% Fix physical mount *.tpoly file
    
    % From 3-year obsolete:
    %
    %     [missing header]
    %     3.258;0;3.258;5.906;0;5.906;0;0;
    %     [Tracker divisions;;;;;;;;;;;;;;;]
    %     Div 0;0;0;1.624;0;1.624;0.976;0;0.976;
    
    trckarea = deprecated_importpolygonsfile(mountfile);
    [trcktlbl,trckinfo] = mountinfo(trckarea,type);
    
    trckspecs.mount_model_id = [PREFIX '_' trcktlbl];
    trckspecs.module_id = module.params.module_id;
    trckspecs.mount_type = type;
    trckspecs.axis_offset = '0;0;0';
    % trckspecs.center_height = NaN;
    trckspecs.comments = {LEGEND; trckinfo;' #center_height: NaN'};
    if contains(type,'1a')
        trckspecs.track_limits = '-45.0; 45.0';
        trckspecs.backtracking = 'on';
    end
    mountfile = fullfile(PATH,[trckspecs.mount_model_id '.tpoly']);
    exportpolygonsfile(trckarea,mountfile,trckspecs,true);
    
%% Fix layout *.mounts file

% #n;x;y;z;[tilt;azimuth;bank]; analysed;shadowthrowing
% 0;2382868.655;5938107.216;0.5;20.0;0.0;0.0;/;/

    layout = readtxtfile(layoutfile,{},'CollectOutput',0);
    assert(numel(layout.data) == 9,'Unexpected layout format');

    layout.headers = {'n;x;y;z;tilt;azimuth;bank;group_id;mount_model'};
    
    Ntr = numel(layout.data{1});
    layout.data{1} = layout.data{1} + 1; % fix zero-based indices
    layout.data{8} = ones(Ntr,1);
    layout.data{9} = ones(Ntr,1);
        
    layout.params.layout_id = sprintf('%s_%dx%s',PREFIX,Ntr,type);
    layout.params.coords_system = 'CUSTOM #WGS84 for decimal lat/lon, UTM XXZ for UTM zone XXZ';
    layout.params.coords_center = '0.000000;0.000000 #LAT;LON, required for CUSTOM projection';
    layout.params.coords_rotation = '0.0 #degrees around zenith, required for CUSTOM projection';
    layout.params.mount_model_idx = ['1; ' trcktlbl];
    
    layout.comments = {LEGEND,[num2str(Ntr) trckinfo]};

    layoutfile = fullfile(PATH,[layout.params.layout_id '.mounts']);
    writetxtfile(layoutfile,layout);

%% Fix array-definition *.arrdef file

%   inverter;combiner a;combiner b;mount;string;modul
%   'a';'c';'b';0;'d';0

%   #file_format: PREPROCESSOR_V2.0
%   #arrdef_id: some_id
%   #layout_id: some_other_id
%   #comments
%   #...
%   #inv;box2;box1;str;mount;pos
%   'a';'b';'c';'d';1;1
%   ...

    arrdef = readtxtfile(arrdeffile,'formatstring','%s %s %s %d %s %d','TreatAsEmpty','/',...
        'EmptyValue',NaN,'CollectOutput',0);
    
    % assert(arrdef.headers,'inverter;combiner a;combiner b;mount;string;modul','Unexpected header');
    arrdef.headers = {'inv;box2;box1;str;mount;pos'};

    arrdef.data = arrdef.data([1 3 2 5 4 6]); % fix column order: inv;box2;box1;str;mount;pos
    arrdef.data{5} = arrdef.data{5} + 1;      % fix 0-based mount & module indices 
    arrdef.data{6} = arrdef.data{6} + 1;
 
    Nm = numel(arrdef.data{1});
    [pcoords,sortidx] = sortrows([arrdef.data{5:6}]);
    [~,ia] = unique(pcoords(:,1),'stable');         % index of first occurence of each mount
    Np = max(int32(diff([ia;Nm+1])));               % modules per mount, based on ia
    if any(pcoords(:,2)-(pcoords(:,1)-1)*Np > Np)
        error('readpparraydef:Npcross','readpparraydef: unexpected module numbering convention.')
    else
        arrdef.data{6}(sortidx) = mod(pcoords(:,2)-1,max(Np))+1;
    end
    
    arrdef.comments = {LEGEND};
    arrdef.params.arrdef_id = PREFIX;
    arrdef.params.layout_id = layout.params.layout_id;
    
    % write (provisional) file
    arrdeffile = fullfile(PATH,[arrdef.params.arrdef_id '.arrdef']);
    writetxtfile(arrdeffile,arrdef);
    lastwill = onCleanup(@() delete(arrdeffile));
    
    % re-read, performing basic checks
    A = readarraydefinition(arrdeffile);
    Ni = max(A(:,3));
    Ns = max(A(:,2));
    Nm = max(A(:,1));
    
    arrdef.params.arrdef_id = sprintf('%s_%di%ds%dm',PREFIX,Ni,Ns,Nm);
    arrdef.comments{end+1} = sprintf('%d inverters with up to %d strings, with up to %d modules',Ni,Ns,Nm);
    arrdeffile = fullfile(PATH,[arrdef.params.arrdef_id '.arrdef']);
    writetxtfile(arrdeffile,arrdef);

end

function [info,modulelbl] = moduleinfo(module)
% Generate module short-label and description
% modulelbl: '980x1960mm_72c_3bd'
% info: '72-cell (3-bypass-diode) 980x1960mm module'

    nc = numel(module.data{1});
    nd = numel(unique(module.data{2}));
    mdims = cellfun(@str2num,strsplit(module.params.outline_polygon,';'));
    mdims = [max(mdims(1:2:end)) - min(mdims(1:2:end)),max(mdims(2:2:end)) - min(mdims(2:2:end))];

    dimlbl = sprintf('%dx%dmm',round(mdims(1)*1000),round(mdims(2)*1000)); % e.g. 980x1960mm
    modulelbl = sprintf('%s_%dc_%dbd',dimlbl,nc,nd);
    
    if nd > 0
        info = sprintf('%d-cell (%d-bypass-diode) %s module',nc,nd,dimlbl);
    else
        info = sprintf('%d-cell %s module (no bypass-diodes)',nc,dimlbl);
    end
end
    
function [mountlbl,info] = mountinfo(A,type)
% Generate mount short-label and description
% mountlbl: '0a_1R20VM_980x1960mm'
% info: 'Fixed tables, 1 row of 20 vertical modules (980x1960mm)'

    P = [A.elements.border];
    mdims = [max(P(1).x) - min(P(1).x),max(P(1).y) - min(P(1).y)];
    if mdims(1) < mdims(2), side = 'vertical'; else, side = 'horizontal'; end
    dimlbl = sprintf('%dx%dmm',round(mdims(1)*1000),round(mdims(2)*1000));   % 980x1960mm

    c = arrayfun(@centroid,P,'unif',0);
    c = cat(2,c{:});
    m = numel(uniquetol(c(1,:),5e-2),'DataScale',1);
    n = numel(uniquetol(c(2,:),5e-2),'DataScale',1);
    if m*n == A.dims
        mountlbl = sprintf('%dR%d%cM',n,m,upper(side(1)));     % 1R20VM
    else
        mountlbl = sprintf('%d%cM',A.dims,upper(side(1)));     % 20VM
    end
    mountlbl = strjoin({type,mountlbl,dimlbl},'_');                          % 0a_1R20VM_980x1960mm
    
    switch type
        case '0a' , typestr = 'Fixed tables';
        case '1aC', typestr = 'Horizontal(contour)-axis trackers';
        case '1aF', typestr = 'Fixed-tilt single-axis trackers';
        case '1aV', typestr = 'Vertical-axis trackers';
        case '2a' , typestr = 'Two-axis trackers';
    end
    
    if m*n == 1
    % e.g. 'Fixed tables, single module (980x1960mm)'
        info = sprintf('%s, single module (%s)',typestr,dimlbl);
    elseif m*n == A.dims
    % e.g. 'Fixed tables, 1 row of 20 vertical modules (980x1960mm)'
        info = sprintf('%s, %s of %s (%s)',typestr,nthings(n,'row'),nthings(m,[side ' module']),dimlbl);
    else
    % e.g. 'Fixed tables (irregular array), 20 vertical modules (980x1960mm)'
        info = sprintf('%s (irregular array), %d %s modules (%s)',typestr,m,side,dimlbl);
    end
end

function writetxtfile(filename,datastruct,format)

    if nargin < 3, format = formatstring(datastruct.data); end

    fileID = fopen(filename,'w');
	assert(fileID >= 0,'exportpolygonsfile: could not create/open file');
    lastwill = onCleanup(@() fclose(fileID));
    
    % Write header
    fprintf(fileID,'#file_format: PREPROCESSOR_V2.0\r\n');
    if isfield(datastruct,'comments')
        cellfun(@(r) fprintf(fileID,'#  %s\r\n',r),datastruct.comments); 
    end
    if isfield(datastruct,'params')
        for f = fieldnames(datastruct.params)'
            fprintf(fileID,'#%s: %s\r\n',f{1},datastruct.params.(f{1})); 
        end
    end
    if isfield(datastruct,'headers')
        cellfun(@(r) fprintf(fileID,'#%s\r\n',r),datastruct.headers); 
    end

    % Write data
    Nt = numel(datastruct.data{1});
    for j = 1:numel(datastruct.data)
       if ~iscell(datastruct.data{j}), datastruct.data{j} = mat2cell(datastruct.data{j},ones(Nt,1)); end 
    end
    datastruct.data = cat(2,datastruct.data{:});
    
    for j = 1:Nt
        fprintf(fileID,[format '\r\n'],datastruct.data{j,:});
    end
end

function PolyAreas = deprecated_importpolygonsfile(polygonsfile)
% Reads a .csv file generated by the SHADE tool containing the definition of a set of polygons and
% their outline. The function can be used to import the tracker outline with module layout, or the
% module outline with cell layout.
% - polygonsfile is a string containing the path and filename of the file to be imported. If
%	ommited, a UI import dialog is open to select a file. The format of the file must be exactly
%	the following:
%
% 		[Outline Polygon;]
% 		Px(0);Py(0);Px(1);Py(1);Px(2);Py(2); ... ;Px(n);Py(n);
% 		Tracker divisions;
% 		Div 0;p1x(0);p1y(0);p1x(1);p1y(1);p1x(2);p1y(2); ... ;p1x(m);p1y(m);
% 		Div 1;p2x(0);p2y(0);p2x(1);p2y(1);p2x(2);p2y(2); ... ;p2x(m);p2y(m);
% 		...
% 		Div k;pkx(0);pky(0);pkx(1);pky(1);pkx(2);pky(2); ... ;pkx(m);pky(m);
% 
% 	Where Px(1..n),Py(1..n) define the n vertexes of the outline polygon, and each row
% 	pix(1..m), piy(1..m) defines m vertexes for the i'th internal polygon.
%
% - PolyAreas is an object of the class pvArea, with two levels depth, containing the outline and
%	subdivision polygons in PolyAreas.border and PolyAreas.elements(j).border, respectively.

	fileID = fopen(polygonsfile);
    assert(fileID > 0,'importpolygonsfile: cannot open selected file');
    
    % Get polygon outline from the 1st or 2nd line
    polygoninfo = textscan(fileID,'%[^\n\r]',1); % read the first line
    hasheader = false;
	try %#ok
        outlinepolygon = readpolygon(polygoninfo{1});
    catch
        polygoninfo = textscan(fileID,'%[^\n\r]',1); % read the 2nd line (outline polygon)
        outlinepolygon = readpolygon(polygoninfo{1});
        hasheader = true;
    end
	
	textscan(fileID,'%*[^\n\r]',1); %skip one line: 'Tracker divisions' ...
	% ... and read the first column, just to count the number of elements
	divisions = textscan(fileID,'Div %d%[^\n\r]');
	divisions = numel(divisions{1});
	
	% initialize a pvArea object with outline polygon and the required number of elements
	PolyAreas = pvArea(outlinepolygon.x,outlinepolygon.y);
	PolyAreas.elements(divisions,1) = pvArea();
	
	frewind(fileID);                            % go back to the beginning...
	textscan(fileID,'%*[^\n\r]',2+hasheader);	%...  and skip the header lines
    for j = 1:divisions                         % then create a pvArea for each polygon
		polygoninfo = textscan(fileID,'Div %d;%[^\n\r]',1);
		thiselement = readpolygon(polygoninfo{2});
		PolyAreas.elements(j) = pvArea(thiselement.x,thiselement.y);
    end
    
    function Poly = readpolygon(polystring)
        % Read a string of the form N;px(0);py(0);px(1);py(1); ... ;px(m);py(m);
        % and interpret it as a polygon with m vertexes px(1..m),py(1..m)

        vertices = textscan(polystring{1},'%f','Delimiter',';');
        vertices = vertices{1};
        if mod(numel(vertices),2) ~= 0, vertices = vertices(2:end); end
        Poly = polygon(vertices(1:2:end),vertices(2:2:end));
    end
end

function s = formatstring(data)
    decimals = @(x) round(8-log10(max(abs(x(isfinite(x))))));

    s = '';
    for j = 1:numel(data)
        v = data{j};
        
        if iscellstr(v) %#ok<ISCLSTR>
            s = [s '%s;']; %#ok<AGROW>
            continue;
        end
        
        if all(mod(v,1) == 0)
        % zero-padded integers
            d = numel(num2str(max(abs(v))));
            s = [s '%0' num2str(d) 'd;']; %#ok<AGROW>
           continue;
        end

        % choose number of decimals ...
        d = decimals(v);
        if ~isfinite(d), d = 2; 
        else
            while all(any(round(mod(v*10^(d-1),1),1) == [0,1],2)), d = d-1; end
        end
        %if abs(d) <= 5
            d = max(0,d);
            s = [s '%0.' num2str(d) 'f;']; %#ok<AGROW>
        % else  
        % % ... or scientific notation
        %     d = decimals(1);
        %     formatstr = [formatstr ';%0.' num2str(d) 'e']; %#ok<AGROW>
        % end
    end
end
		