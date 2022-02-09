function ModuleArea = stdmodulepoly(varargin)
% Generates a pvarea object for a rectangular module with square cells.
%
% STDMODULEPOLY() - 72 cell (3 diode) horizontal module with default size and margin
% STDMODULEPOLY([nx,ny]) - nx·ny cell module with default size and margins
%                             if mod(min(nx,ny),3)=0, module is split lengthwise in 3 substrings.
% STDMODULEPOLY([nx,ny,ndx,ndy]) - module is split in ndx·ndy substrings.
% STDMODULEPOLY([nx,ny,...],[w,h]) - nx·ny cell module with margins adjusted to march size w·h (*)
% STDMODULEPOLY([nx,ny,...],[w,h],gap) - just as above, with custom gap [gx,gy] between cells  (*)
% STDMODULEPOLY([nx,ny,...],[l,r,t,b]) - nx·ny cell module with margins [l,r,t,b]
% STDMODULEPOLY([nx,ny,...],[l,r,t,b],gap) - ...
% STDMODULEPOLY([nx,ny,...],[l,r,t,b],gap,cellsize) - ... you get the idea
%
% STDMODULEPOLY([],[w,h]) - try to find a 156 mm-cell, 6 - row module that fits the provided   (*)
%                           dimensions with reasonable gaps and margins.
% STDMODULEPOLY(obj) - where obj is a pvArea, polygon object, or ODM structure of parameters will
%   try to do the same.        
%
% (*) NOTE: module dimensions will always be rounded to 2·mm, so that centered border coordinates
%           can be written with three decimal places.
%
% stdmodulepoly(...,'outputfile',filename) generates file filename.mpoly in cd.
%
% [nx,ny] - number of cells in horizontal and vertical direction
% [..,ndx,ndy] - number of substrings (diodes) along each direction
% [w,h] - module width and height (defaults calculated from gaps and margins, see below)
% [l,r,t,b] - module margins
% gap [gx,gy] - gaps between cells
% cellsize [cx,cy]
%
% See also: EXPORTPOLYGONSFILE, SAMPLETRACKERS

%  DEFAULTS (vertical vectors!):

    OPT.ncells = [12;6];                 % Horizontal 72 cell module
    OPT.cellsize = [1;1]*0.156;          % 156 mm square cells
    OPT.margins = [1;1;1;1]*0.12*0.156;  % 12% of cell size in all directions
    OPT.gap = [1;1]*0.03*0.156;          % 3% of cell size in both directions

	OPT.maxgap = OPT.gap*1.5;            % limits used to guess number of cells
    OPT.mingap = OPT.gap/1.5; 
    OPT.maxmarg = OPT.margins*1.5;
    OPT.minmarg = OPT.margins/1.5;
    
    OPT.outputfile = '';
    OPT.ok2write = [];
    
    %Look for outputfile option
    [OPT, args] = getpairedoptions(varargin,OPT);
    
    % Remove trailing empty-arguments
    while ~isempty(args) && isempty(args{end}), args = args(1:end-1); end

    switch numel(args)
        case 0     
        % stdmodulepoly() - use default cells and margins
            ModuleArea = stdmodulepoly(OPT.ncells,OPT.margins);
        case 1     
            % stdmodulepoly(obj) - switch to stdmodulepoly([],[w,h])
            switch class(args{1})
                case 'struct'
                    S = args{1};
                    assert(all(isfield(S,{'Ns','width','length'})),...
                                            'Expected fields {Ns, width, length}');
                    if isfield(S,'Np')
                        assert(S.Np == 1,'Cannot yet handle modules with Np ~= 1');
                    end
                    
                    dims = sort([S.width,S.length]);
                    if prod(dims) > 1e5, dims = dims/1000;
                    elseif prod(dims) > 1e3, dims = dims/100;
                    end
                    ncells = guessncells(OPT,dims);
                    assert(S.Ns == prod(ncells(1:2)),'Unknown cell configuration');
                    if isfield(S,'Nbpd')
                    % Try dividing the shortest side... out of ideas here
                        [~,short] = min(ncells);
                        assert(mod(ncells(short),S.Nbpd) == 0,'Unknown bypass configuration');
                        ncells(3:4) = 1;
                        ncells(find(short)+2) = S.Nbpd;
                    end
                    ModuleArea = stdmodulepoly(ncells,dims);
                case 'polygon'
                    [w,h] = rectangleproperties(args{1});
                    ModuleArea = stdmodulepoly([],[w,h]);
                case 'pvArea'
                    [w,h] = rectangleproperties(args{1}.border);
                    ModuleArea = stdmodulepoly([],[w,h]);
                otherwise
                % stdmodulepoly([nx,ny,...])
                    assert(mod(varargin{1},1) == 0,'expecting integer number of cells');
                    ModuleArea = stdmodulepoly(varargin{1}(:),OPT.margins); % set default margins and move on
            end
        case 2
            switch numel(varargin{2})
                case 4  
                % stdmodulepoly([nx,ny,...],[l,r,t,b])
                    % set default cell gap and move on
                    ModuleArea = stdmodulepoly(varargin{1}(:),varargin{2}(:),OPT.gap);
                    
                case 2  
                % stdmodulepoly([...],[w,h])

                    % see if we know the number of cells, try to guess otherwise
                    if isempty(varargin{1}) 
                    % stdmodulepoly([],[w,h])
                       ncells = guessncells(OPT,varargin{2}); 
                    else
                        ncells = varargin{1}(:);
                    end
                    
                    % stdmodulepoly([nx,ny,...],[w,h])
                    ncells = ncells(1:2);       % ignore substring distribution for now**
                    cellsize = OPT.cellsize;
                    moddims = varargin{2}(:);
                    
                    assert(prod(moddims) < 100,'You reeally want to check your module-dimensions');
                    
                    % check if the default cells will fit in the provided dimensions,
                    % and issue warnings if gaps or margins fall outside the DEF thresholds
                    space = moddims - ncells.*cellsize;

                    if any (ncells == 0) || any(space < 0)
                        error('stdmodulepoly:cells do not seem to fit in module'); 
                    end

                    if any((ncells-1).*OPT.mingap + sum(reshape(OPT.minmarg,2,2))' > space)
                        warning('stdmodulepoly: cell gaps and/or module margins are narrower than normal')
                    elseif any((ncells-1).*OPT.maxgap + sum(reshape(OPT.maxmarg,2,2))' < space)
                        warning('stdmodulepoly: cell gaps and/or module margins are larger than normal')
                    end

                    % divide the available space proportionally to default gaps and margins
                    factor = space./((ncells-1).*OPT.gap + sum(reshape(OPT.margins,2,2))');
                    
                    gap = OPT.gap.*factor;
                    margins = OPT.margins.*reshape([factor,factor]',4,1);
                    
                    % Put substring information back on, in case we had removed it**
                    if numel(varargin{1})> 2, ncells = varargin{1}; end
                    
                    ModuleArea = stdmodulepoly(ncells,margins,gap,cellsize);
                otherwise
                    error('stdmodulepoly:cannot regognize arguments');
            end
        case 3  % stdmodulepoly([...],[...],gap
            gap = varargin{3}(:); if isscalar(gap), gap = [gap;gap]; end
            
            if any(gap < OPT.mingap)
                warning('stdmodulepoly: cell gaps are narrower than normal')
            elseif any(gap > OPT.maxgap)
                warning('stdmodulepoly: cell gaps are larger than normal')
            end
            
            switch numel(varargin{2})
                case 2
                % stdmodulepoly([...],[w,h],gap
                
                    % see if we know the number of cells, try to guess otherwise
                    if isempty(varargin{1}) 
                    % stdmodulepoly([],[w,h],gap)
                       ncells = guessncells(OPT,varargin{2},gap); 
                    else
                        ncells = varargin{1}(:);
                    end
                    
                    % stdmodulepoly([nx,ny,...],[w,h],gap)
                    
                    ncells = ncells(1:2);       % ignore substring distribution for now**
                    cellsize = OPT.cellsize;
                    moddims = varargin{2}(:);
                    
                    % check if the default cells will fit in the provided dimensions,
                    % and issue warnings if gaps or margins fall outside the DEF thresholds
                    space = moddims - ncells.*cellsize-(ncells-1).*gap;

                    if any(space < 0)
                        error('stdmodulepoly:cells do not seem to fit in module using the provided gap'); 
                    end

                    if any(sum(reshape(OPT.minmarg,2,2))' > space)
                        warning('stdmodulepoly: module margins are narrower than normal')
                    elseif any(sum(reshape(OPT.maxmarg,2,2))' < space)
                        warning('stdmodulepoly: module margins are larger than normal')
                    end
    
                    % Center the cell array
                    margins = reshape([space,space]',4,1)/2;
                    
                    % Put substring information back on, in case we had removed it**
                    if numel(varargin{1})> 2, ncells = varargin{1}; end
                    
                    ModuleArea = stdmodulepoly(ncells,margins,gap,cellsize);
                
                case 4
                % stdmodulepoly([nx,ny,...],[l,r,t,b],gap
                    
                    ModuleArea = stdmodulepoly(varargin{1}(:),varargin{2}(:),gap,OPT.cellsize);
                otherwise
                    error('stdmodulepoly:cannot regognize arguments');
            end
        case 4  
        % stdmodulepoly([nx,ny,...],[l,r,t,b],gap,cellsize) - do the actual work
       
            ncells = varargin{1}(:);
            
            if numel(ncells)>2
               diodes = ncells(3:end);
               ncells = ncells(1:2);
               if any(diodes<1) || any(mod(ncells./diodes,1)) % if division don't make sense
                   error('stdmodulepoly:cannot divide cells in requested substrings');
               end
            else
              if mod(min(ncells),3)==0                      % if short side is divisible by 3
                 diodes = [1;1]+(ncells==min(ncells))*2;    % split lengthwise in 3 strings
              else
                 diodes = [1;1];                            % otherwise use a single substring
              end
            end
            
            margins = varargin{2}(:);
            gap = varargin{3}(:);
            cellsize = varargin{4}(:);
            margins = reshape(margins,[2,2])';     %[l r; t b] - so that all xx and yy share one line
            
            % (*) Round module-dimensions to 2-mm (let margins absorb the difference)
            moddims = ncells.*cellsize + (ncells-1).*gap + sum(margins,2); 
            margins = bsxfun(@plus,margins,(round(moddims/2,3) - moddims/2));
                        
            CellArea = pvArea(cellsize(1),cellsize(2));                              % make single-cell pvArea                    
            stringdims = ncells./diodes.*cellsize + (ncells./diodes-1).*gap;         % substring size
            step = gap + cellsize; 
            cx = cellsize(1)/2+(0:ncells(1)/diodes(1)-1)*step(1)-stringdims(1)/2;    % cell-center vectors for a substring
            cy = cellsize(2)/2+(0:ncells(2)/diodes(2)-1)*step(2)-stringdims(2)/2;

            % arrange cell numbers along the long side
            if ncells(1)/diodes(1) < ncells(2)/diodes(2)
                [cx,cy] = meshgrid(cx,cy);
            else
                [cy,cx] = meshgrid(cy,cx);
            end                        
            CellString = areaArray(CellArea,cx,cy);    % Create the substring as an array of cells
            
            moddims = ncells.*cellsize + (ncells-1).*gap + sum(margins,2);           % module size    
            step = gap + stringdims; 
            cx = margins(1)+stringdims(1)/2+(0:diodes(1)-1)*step(1)-moddims(1)/2;    % string-center vectors
            cy = margins(4)+stringdims(2)/2+(0:diodes(2)-1)*step(2)-moddims(2)/2;
            
            % arrange string numbers along the long side
            if diodes(1) < diodes(2),[cx,cy] = meshgrid(cx,cy);
            else, [cy,cx] = meshgrid(cy,cx); end
            
            ModuleArea = areaArray(CellString,cx,cy);    % Create the module as an array of cell-strings
            ModuleArea.border = polygon(moddims(1),moddims(2));    % Increase border area to consider margins
            
            %figure; plotArea(ModuleArea,[],true);
            
        otherwise
             error('stdmodulepoly:don''t know what to do with so many arguments');
    end
    
    % Generate module files
    if ~isempty(OPT.outputfile)
        specs.module_id = regexprep(OPT.outputfile,'\.[^.]{0,10}$',''); % remove extension   
        specs.comments = {'Generated with STDMODULEPOLY'};
        
        filename = [specs.module_id '.mpoly'];

        % Check if files exist, and if so, if they can be overwritten
        if isempty(OPT.ok2write), OPT.ok2write = right2overwrite(filename); end
        if ~OPT.ok2write; return; end
        
        exportpolygonsfile(ModuleArea,filename,specs,OPT.ok2write);        % write module polygons
    end
end

function ncells = guessncells(DEF,moddims,gap)
% Try to guess the number of cells and cell-strings [nx ny ndx ndy] that fit in a module of 
% dimensions dims [w,h].

    % Recognized 'standard' cell numbers:
    stdsizes = [12,8; 12,6; 10,6; 9,6; 9,4; 8,4; 6,6]';
    
    if nargin < 3, gap = DEF.gap; end
    
    ncells = (moddims(:)+ gap - sum(reshape(DEF.margins,2,2))')./(DEF.cellsize + gap);
    ncells = round(ncells);
    
    stdsizes = [stdsizes,flipud(stdsizes)]; % take also rotated module
    
    if ~any(all(bsxfun(@eq,ncells,stdsizes),1))
       warning('stdmodulepoly:guessncells','%d x %d cells is not a recognized module standard', ncells(1),ncells(2));
    end
end
