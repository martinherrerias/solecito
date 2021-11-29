function [MergeGeom,flag] = fitmoduleinmount(ModuleGeom,MountGeom)
% Replace each ModuleGeom.elements(j) with a translated (and rotated) copy of MountGeom
% check for compatibility and overlap issues.

    TOL = 0.05; % should be larger than equal-element tolerance in IMPORTPOLYGONSFILE

    % Quick-check Module area against Tracker-element(s) area
    areaoffset = [MountGeom.elements.area]/ModuleGeom.area - 1;
    if any(abs(areaoffset) > TOL)
        switch any(areaoffset > 0) - any(areaoffset < 0)
            case 1, sgn = 'larger';
            case 0, sgn = 'larger/smaller';
            case -1, sgn = 'smaller';
        end
        warning('GUIphysmod:mountmodfit','Mount-element-area(s) are up to %0.1f%% %s than the provided module''s area',...
            max(abs(areaoffset)),sgn);
    end

    % Remove any offset from centroid
    ModuleGeom = translate(ModuleGeom,-centroid(ModuleGeom.border));
    
    % Detect 90° rotation of modules
    % FUTURE: instead of providing outline polygons, MOUNT-FILES should just be a list of centroids
    % and rotation angles, allowing non-symetric modules (connection pads) and ~90° angles.
    modflip = arrayfun(@moduleorientation,[MountGeom.elements.border]);
    defflip = moduleorientation(ModuleGeom.border);
    if any(modflip == 'u') || defflip == 'u'
        warning('Double-check the orientation of the module-definition, as we can''t make it out from the mount-polygons')
        modflip(modflip == 'u') = 'h';
        if defflip == 'u', defflip = 'h'; end
    end
    rot = 90*(modflip ~= defflip)';
    
    c = cell2mat(arrayfun(@centroid,[MountGeom.elements.border],'unif',0))';
    
    % Make sure ModuleAreas is 3-levels-deep (module-group-cell)
    % so that ModuleAreas.element(1).elements(1) returns a cell
    if ModuleGeom.depth ~= 3
        error('pvArea object of depth 3 required')
        % for j = 1:ModuleGeom.size
        %    ModuleGeom.elements(j).elements(1) = copy(ModuleGeom.elements(j)); 
        % end
    end
    
    % Replicate modules onto tracker
    MergeGeom = pvArea(MountGeom.border,...
        arrayfun(@(a,x,y) translate(rotate(copy(ModuleGeom),a),x,y),rot,c(:,1),c(:,2)));
    checkoverlaps(MergeGeom);
    flag = 0;
end

function key = moduleorientation(p)
    [w,h,~] = rectangleproperties(p);
    if w == h, key = 'u';       % unknown
    elseif w > h, key = 'h';    % horizontal
    else, key = 'v';            % vertical
    end
end

function checkoverlaps(obj)
% Check if any two elements of pvArea object OBJ intersect at some point.
% Issue an error if they do.

    P = [obj.elements.border];
    N = numel(P);
    if N < 2, return; end

    P = pack(P); % pre-pack for speed and numerical-consistency
    for k = N:-1:1
       q = polyclip(P(k),P([1:k-1,k+1:end]),1,'pos','pos'); % intersection of k with all-but-k
       overlap(k) = polygon.packedarea(q)/polygon.packedarea(P(k));
    end
    
    if any(overlap > 0)
        error('GUIphysmod:checkoverlaps','%s have overlaps of up to %0.2f%%',...
            shortliststr(find(overlap > 0),'Element',10),max(overlap));
    end
end