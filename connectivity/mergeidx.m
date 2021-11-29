function C = mergeidx(B,A,mergeorphans)
% C = MERGEIDX(B,A) 
% C = MERGEIDX(B,A,mergeorphans)

% MERGEIDX is designed to merge two NDmap/LabelMap-structures with fields representing maps in
%   one (or several) dimension(s). The resulting structure C will have:
%   - All the common fields of A and B, for which C.x = mergemap(B.x,A.x). 
%   - All 'composite' fields (those which can be formed by fields in the other structure) 
%     e.g. C.xyz = mergemap(compose(B.xy,B.z),A.xyz), this is useful when some dimensions are
%     'linked' in one index, but not in the other.
%   - If mergeorphans = true (default is false), any NDmaps from either A or B which are not
%     included in the points above (also not as part of a composite) will just be copied to C.
%
% See also: MERGEMAP, NDMAP, LABELMAP

    if nargin < 3, mergeorphans = false; end

    fa = fieldnames(A);
    fb = fieldnames(B);
    
    % Merge common fields of A and B
    f = intersect(fa,fb);
    for j = 1:numel(f)
        C.(f{j}) = mergemap(B.(f{j}),A.(f{j}));
    end
    
    % Try to merge 'composite' fields from both A and B
    mergedfb = mergecomposites(A,B,fa,fb,false);
    mergedfa = mergecomposites(B,A,fb,fa,true);
    
    % Copy 'orphan' fields from A and B
    if mergeorphans
        fc = fieldnames(C);
        f = setdiff(setdiff(fa,mergedfa),fc);
        for j = 1:numel(f), C.(f{j}) = A.(f{j}); end
        f = setdiff(setdiff(fb,mergedfb),fc);
        for j = 1:numel(f), C.(f{j}) = B.(f{j}); end
    end
    
    function usedfb = mergecomposites(A,B,fa,fb,flipped)
    % Try to find fields in A which can be composed from 'unit' fields in B
    
        usedfb = {};
        fa = setdiff(fa,fb);
        for k = 1:numel(fa)
            if numel(fa{k}) <= 1, continue; end % skip unitary fields

            p = findpieces(fa{k},fb); % try to get 'unit' fields in B that compose fk in A
            if isempty(p), continue; end
            
            if iscell(p{1})
                warning('mergeidx:compose',...
                    'field ''%s'' can be composed in more than one way, using: %s',...
                    fa{k}, strjoin(cellfun(@(x) sprintf('''%s''',x),p{1},'unif',0),', ')); 
                p = p{1};
            end

            % Dump 'unitary' B-structure fields p into a cell array of arguments
            args = cellfun(@(x) B.(x),p,'unif',0);
            usedfb = [usedfb{:},p];

            % Create a composite NDmap from 'unit' fields in B, and merge with A.(fa{k})
            if ~flipped
                C.(fa{k}) = mergemap(compose(args{:}),A.(fa{k}));
            else
                C.(fa{k}) = mergemap(A.(fa{k}),compose(args{:}));
            end
        end
    end
end

function p = findpieces(str,pieces)
% p = FINDPIECES(str,pieces) Quick-and-dirty approach to find combination(s) of pieces 
% (cell-array-of-strings) that, put together, form a string str. e.g. 
%
%   p = findpieces('abcd',{'a','d','x','bc'})       -> p = {'a','bc','d'}
%   p = findpieces('abcd',{'a','d','x','bc','abc'}) -> p = {{'a','bc','d'},{'abc','d'}}
%   p = findpieces('nope',{'a','b','c'})            -> p = {}
%
% CAUTION: the function evaluates all combinations of pieces which appear somewhere in str, and 
%   all permutations of combinations with propper length. i.e. things can get quickly out-of-hand.
%   Use only for short strings and a handfull of pieces.

    p = {};

    % Stick to pieces that appear somewhere in str
    useful = pieces(cellfun(@(c) contains(str,c),pieces));
    if isempty(useful), return; end

    s = cellfun(@numel,useful); % vector of piece-sizes

    % Get min and max number of pieces required to form str
    k(2) = min(numel(str),ceil(numel(str)/min(s)));
    k(1) = max(1,floor(numel(str)/max(s)));

    % Get a list of indices for all possible combinations with k(1)-to-k(2) useful-pieces 
    c = arrayfun(@(j) num2cell(nchoosek(1:numel(useful),j),2),k(1):k(2),'unif',0); 
    c = cat(1,c{:});

    % Keep only those whith the correct size
    rightsize = cellfun(@(x) sum(s(x)),c) == numel(str);
    c = c(rightsize);
    if isempty(c), return; end

    % Consider all possible permutations
    c = cellfun(@(x) num2cell(perms(x),2),c,'unif',0);
    c = cat(1,c{:});

    % Get the actual composed strings, and keep those that work
    cc = cellfun(@(x) cat(2,useful{x}),c,'unif',0);
    c = c(strcmp(str,cc));

    switch numel(c)
        case 0
        case 1
            p = useful(c{1})';
        otherwise
            p = cellfun(@(x) useful(x)',c,'unif',0);
    end
end
