function exportpolygonsfile(A,filename,specs,forceoverwrite)
% EXPORTPOLYGONSFILE(A,FILENAME,SPECS)
% EXPORTPOLYGONSFILE(...,FORCEOVERWRITE)
% Creates a *.mpoly/*.tpoly file describing a module- or mount(tracker)-geometry, including any 
% specifications passed in the SPECS structure, and a (default) interconnection scheme extracted 
% from the structure of PVAREA A. Output format follows the 2.0 Preprocessor Specification:
%
%     #file_format: PREPROCESSOR_V2.0
%     #SPECS_field_1: value_1
%     #SPECS_field_2: value_2
%     #...
%     #outline_polygon: Px(1);Py(1);Px(2);Py(2); ... ;Px(n);Py(n);
%     #divisions(n;string_id;x1;y1;x2;y2;...)
%     1;1;p1x(1);p1y(1);p1x(2);p1y(2); ... ;p1x(m);p1y(m);
%     2;1;p2x(1);p2y(1);p2x(2);p2y(2); ... ;p2x(m);p2y(m);
%     ...
%     q;1;pqx(1);pqy(1);pqx(2);pqy(2); ... ;pqx(m);pqy(m);
%
% 	Where Px(1..n),Py(1..n) define the n vertices of the outline polygon, and each row
% 	i,j,pix(1..m), piy(1..m) defines the i'th subdivision (cell/module), as an m-vertex polygon
%   connected to (cell-) string j.
%
% A: is an object of the class pvArea, with two or three levels depth, containing the outline 
%   (A.border) and subdivision polygons A[.elements(k)].elements(j).border. For depth-3 objects,
%   the intermediate level A.elements(k) groups subdivisions into strings.
%
% FILENAME: string, including extension
%
% For module-model definitions, SPECS should include the following fields:
%  SPECS.module_id: (string) module ID
%  SPECS.comments: optional cell-array of strings
% 
% For mount-model definitions, depending on mount-type SPECS might include the following fields:
%  SPECS.mount_model_id: (string) Mount-Model ID
%  SPECS.module_id: (string) Module ID that matches an existing module_id.mpoly file
%  SPECS.mount_type: (string) {2a,1aV,1aF,1aC,0a}
%  SPECS.axis_offset: 3-vector
%  SPECS.center_height: scalar
%  SPECS.track_limits: 2/4-vector, depending on tracker type
%  SPECS.backtracking: (string) {on,off,smart}
%  SPECS.comments: optional cell-array of strings
%
% FORCEOVERWRITE: boolean, overwrite any existing FILENAME
%
% See also: WRITETRACKERLAYOUT, EXPORTPOLYGONSFILE, WRITEARRAYDEFINITION

    narginchk(1,4);
	if nargin < 2 || isempty(filename), filename = 'export_poly.txt'; end
    if nargin < 3 || isempty(specs), specs = struct('comments',{'Unknown object type'}); end
    if nargin < 4, forceoverwrite = false; end
		
    if ~right2overwrite(filename,forceoverwrite), return; end
    assert(isa(A,'pvArea') && A.depth > 1 && A.depth < 4,'A must be a 2 to 3-depth pvArea');
    
 	fileID = fopen(filename,'w');
	assert(fileID >= 0,'exportpolygonsfile: could not create/open file');

    % Write header
    fprintf(fileID,'#file_format: PREPROCESSOR_V2.0\r\n');
    specfields = setdiff(fieldnames(specs),{'comments'},'stable');
    for f = specfields'
        fprintf(fileID,'#%s: %s\r\n',f{1},specs.(f{1})); 
    end
    if isfield(specs,'comments')
        cellfun(@(r) fprintf(fileID,'#  %s\r\n',r),specs.comments); 
    end
    
    % Write outline
	fprintf(fileID,'#outline_polygon: %s\r\n',writepolygon(A.border));
    
	% Write divisions
	fprintf(fileID,'#divisions(n;string_id;x1;y1;x2;y2;...)\r\n');
    if A.depth == 2
        for j = 1:numel(A.elements)
            fprintf(fileID,'%d;1;%s\r\n',j,writepolygon(A.elements(j).border));
        end
    else
        idx = 1;
        for k = 1:numel(A.elements)
            for j = 1:numel(A.elements(k).elements)
                fprintf(fileID,'%d;%d;%s\r\n',idx,k,writepolygon(A.elements(k).elements(j).border));
                idx = idx + 1;
            end
        end
    end
	fclose(fileID);
end

function polystring = writepolygon(Poly)
% Print the set of m-vertices of polygon POLY as a string of the form:
% pqx(1);pqy(1);pqx(2);pqy(2); ... ;pqx(m);pqy(m);
	polystring = sprintf('%0.4f;%0.4f;',[Poly.x(:),Poly.y(:)]');
end
		
