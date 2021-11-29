function [xc,yc] = proj(prj4,x,y)
% [XC,YC] = PROJ(PRJ4,X,Y) - wrap for OSGeo/proj coordinate transformation software library [1]
%
% [xp,yp] = PROJ('+proj=...',lon,lat) - convert [lat,lon] using PRJ.4 projection '+proj...'
% [lon,lat] = PROJ('-I +proj=...',xp,yp) - inverse PRJ.4 projection '+proj...'
% [xc,yc] = PROJ(C,X,Y) - for cellstr C of size(X), uses a different projection for each* point,
%       (*) Internally, data is clustered by unique(C), to reduce OSGeo/proj calls.
%
%   Examples:
%   [x,y] = proj('+proj=utm +zone=32',9.097,48.74);
%   [a,b] = proj('-I +proj=tmerc +lat_0=48.74 +lon_0=9.097 +x_0=0 +y_0=0',100,50);
%   proj() % run proj.test, an interactive ISEA world projection example
%
% [1] PROJ contributors (2020). PROJ coordinate transformation software library. 
% Open Source Geospatial Foundation. URL https://proj.org/.
%
% See also: UTM2DEG, DEG2UTM, PRJ2ABS, ABS2PRJ

    if nargin == 0, test(); return; end
    narginchk(3,3);
    validateattributes(x,{'numeric'},{'real'},'proj','x');
    validateattributes(y,{'numeric'},{'real','size',size(x)},'proj','y');

    if iscellstr(prj4) || isstring(prj4)
        if isscalar(prj4)
            prj4 = prj4{1};
        else
            assert(numel(prj4) == numel(x),'Expecting scalar or numel(X) PRJ4 string');
            [xc,yc] = batchproj(prj4,x,y);
            return;
        end
    end
    
    assert(ischar(prj4) && contains(prj4,'+proj='),'Expecting char-string +proj=...');
    prj4 = strtrim(prj4);

    % Check for PROJ <https://proj.org/index.html> on system
    [out,~] = system('echo 0 0 | proj +proj=ortho');
    assert(out == 0,'Error testing proj on system');
    
    [out,ARG_MAX] = system('getconf ARG_MAX');
    if out, ARG_MAX = 128*1024; else, ARG_MAX = str2double(ARG_MAX)/16; end

    if contains(prj4,'-I')
        cmd = sprintf('proj -f "%%.8f" -e "NaN\tNaN" ');
        fmt = '%0.8e %0.8e\n';
    else
        cmd = sprintf('proj -f "%%.3f" -e "NaN\tNaN" ');
        fmt = '%0.8f %0.8f\n';
    end
    
    sz = size(x);
    if ~isvector(x) || size(x,2) > 1, x = x(:); y = y(:); end
    
    xc = NaN(sz);
    yc = NaN(sz);   
    ok = isfinite(x) & isfinite(y);
    if ~any(ok), return; end
    
    % Function calls are limited by console characters
    % TODO: proj is supposed to take binary data directly
    chunks = ceil((nnz(ok)*28 + numel(prj4) + 16)/ARG_MAX);
    
%     if chunks == 1
%         file = sprintf(fmt,[x(ok),y(ok)]');
%         [~,out] = system([cmd prj4 ' <<EOF' newline() file 'EOF']);
%         C = str2num(out); %#ok<ST2NM>
%         xc(ok) = C(:,1);
%         yc(ok) = C(:,2);
%         return;
%     end

    x = x(ok);
    y = y(ok);
    C = cell(numel(chunks));
    n = ceil(numel(x)/chunks);    
    inj = 0;
    for j = 1:chunks
        c = inj+1:min(numel(x),inj+n);
        inj = inj+n;
        file = sprintf(fmt,[x(c),y(c)]');
        [~,out] = system([cmd prj4 ' <<EOF' newline() file 'EOF']);
        C{j} = str2num(out); %#ok<ST2NM>
        if isempty(C{j})
            out = strtrim(regexprep(out,'.*<proj>:(.*)','$1'));
            error('OSGeo/proj error: %s',out);
        end
    end
    C = cat(1,C{:});

    xc(ok) = C(:,1);
    yc(ok) = C(:,2);
end

function [xc,yc] = batchproj(prj4,x,y)

    [~,ic,ia] = uniquecell(prj4);
    u = numel(ic);

    xc = NaN(size(x));
    yc = NaN(size(y));

    failed = true(u,1);
    errs = cell(u,1);
    for k = 1:u
        j = ic(k);
        inj = ia == k;
        try
            [xc(inj),yc(inj)] = proj(prj4{j},x(inj),y(inj));
            failed(k) = false;
        catch ERR
            errs{k} = getReport(ERR);
        end
    end
    if any(failed)
        warning('%d of %d projections resulted on %s',nnz(failed),u,...
            shortliststr(errs(failed),'error','-newline'));
    end
end

function test()
% TEST function for both PLOTCONTROLS and PROJ. Draws an icosahedral world projection ('isea')
% with interactive controls (one menu of presets + three sliders).

    % Initialize plot
    load('coastlines','coastlat','coastlon');
    [xp,yp] = proj('+proj=isea +a=0.95492966 +orient=isea',coastlon,coastlat);
    
    GUIfigure('projtest','proj test','3:1'); clf();
    set(subplot(1,4,4),'visible',false);
    subplot(1,4,1:3); 
    hold on;
    set(gca,'visible',false)
    icos = [(-2.5:0.5:2.5)*2/sqrt(3);0.5+[0 1 0 1 0 1 0 1 0 1 0]];
    icos = [icos,icos(:,1),icos + [1/sqrt(3);-1],icos.*[1;-1]+[1/sqrt(3);0]];
    plot(icos(1,:),icos(2,:));
    
    h = scatter(xp,yp,1);
    axis equal

    % Define presets
    PRESETS.isea = {58.28 5.62 0}; % fitted, doesn't seem to match PROJ code (?)
    PRESETS.pole = {0,90,0};
    PRESETS.fuller = {2.3,-5.2454,7.46658}; % found as comment in PROJ code, no idea if true

    % add one menu with presets, and three sliders
    ctrls = plotcontrols({'m','s','s','s'},{'orient','lat_0','lon_0','azi'},...
        {[fieldnames(PRESETS);{'custom'}]',[-90,90,1,5],[-180,180,0.9,9],[-180,180,1,5]},...
        ['isea',PRESETS.isea],@updatewithpreset,'position',subplot(1,4,4).Position,...
        '-skipupdate','-smooth');
    
    disabled = false; % (§)
    
    function plotupdate(varargin)
    % plotupdate(lat_0,lon_0,azi) - reproject and update data points
        prj4 = sprintf('+proj=isea +a=0.95492966 +lat_0=%0.3g +lon_0=%0.3g +azi=%0.3g',...
            varargin{:});   
        [h.XData,h.YData] = proj(prj4,coastlon,coastlat);
        drawnow();
    end

    function updatewithpreset(varargin)
    % Callback for PLOTCONTROLS, called as UPDATEWITHPRESETS(P,A,B,..,SRC,EVENT) where A,B.. 
    % are the values of the actual UPDATE arguments (in this case: lat_0, lon_0, azi)
    %
    % P (varargin{1}) is taken as a preset configuration string set by ctrls(1), which, if set 
    % by the user, determines the values of A,B,..
    %
    % If the update was triggered by another control, ctrls(1) is set to 'custom'
    % If ctrls(1) is set by the user to 'custom', nothing changes
    %
    % TODO: take dynamic lists? with the possibility for the user to 'save' presets?

        if disabled, return; end % (§)

        if isequal(varargin{end}.AffectedObject,ctrls(1))    
            if isfield(PRESETS,varargin{1})
             % Set sliders to preset upon menu change
                varargin(2:end-2) = PRESETS.(varargin{1});
                resetcontrols()
            else
                assert(strcmp(varargin{1},'custom'),'Undefined PRESET?');
            end
        else
        % Set menu to custom upon slider change
            varargin{1} = 'custom'; resetcontrols();
        end

        % Update plot
        plotupdate(varargin{2:end-2});
        
        function resetcontrols()
        % (§) Since ctrls(j).Value = X; will trigger a new call to UPDATE, we have to set a
        % persistent disable while doing any changes to other controls, then return to the update
        
            disabled = true;
            
            % Adjust menu value, in case string (varargin{1}) has been hard set
            [~,ctrls(1).Value] = ismember(varargin{1},ctrls(1).String);
            
            % Adjust slider values, in case they've been set by preset
            for j = 2:numel(ctrls)
               ctrls(j).Value = varargin{j};
            end
            
            disabled = false;
        end
    end
end
