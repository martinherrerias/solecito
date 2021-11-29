function [arrdef,header] = samplearraydef(n,m,N,M,sn,sm,sN,sM,outputfile)
% Generate a simple array definition for an array of N rows of M mounts, each n rows of m modules.
% Modules are connected in strings of sn rows with sm columns, and each sN·sM array of strings
% connects to one inverter/MPPT-input.
% The result represents an [n·N/sN·sn, m·M/sM·sm] array of inverters. Where both n·N/sN·sn
% and m·M/sM·sm must be integers.
%
% arrdef: a 5-column table (mod,str,inv,pos,trck) as per (deprecated) ArrayDefinition convention
% if outputfile is a filename string, the array definition will be written to either
%         a simple CSV file, or a multi-sheet .xls file for readability.
%
% IMPORTANT: it is assumed that physical module layout responds to [yy,xx] = meshgrid(y,x);
% that is, both inside the trackers and in the layout, numbering runs left -> right, down -> up
    
    if nargin < 7 || isempty(sN), sN = N*n/sn; end  % single-inverter as default
    if nargin < 8 || isempty(sM), sM = M*m/sm; end  % ..
    createfile = nargin > 8;

    % Get electrical definition
    br = n*N/(sN*sn);  % inverter-block rows
    bc = m*M/(sM*sm);  % inverter-block columns
    
    isint = @(x) mod(x,1)==0;
    
    if ~isint(br) || ~isint(bc)
        error('samplearraydef:div','Modules not divisible by string blocks')
    end
                
    block = singleblock(sn,sm,sN,sM);   % 2-layer arrdef of single block
    arrdef = repmat(block,br,bc);       % 2-layer arrdef of array
    
    % Get inverter indices:
    arrdef(:,:,[4,3]) = singleblock(sN*sn,sM*sm,br,bc);  % indices at 4th layer will be dumped ..(*)
   
    % Get mechanical definition
    arrdef(:,:,4:5) = singleblock(n,m,N,M); %(*).. here

    % Write array definition to file (before flattening array, in case multi-page is required
    header = {'mod','str','inv','@','pos','trck'};
    if createfile, writearraydefinition(arrdef,outputfile,header); end
    
    % Flatten array definition
    arrdef = shiftdim(arrdef,2);            % shift to a 5·Nn·Mm array
    arrdef = reshape(arrdef,[5,N*n*M*m])';  % turn into a table, N·n·M·m rows by 5 colums
    
    arrdef = sortrows(arrdef,[3,2,1,5,4]);  % sort by inverter > string > module > mount > position
end

function block = singleblock(n,m,N,M)
% Create an [N·n,M·m,2] array of indices, whose first plane is an N·M repetition of n·m 'element'
% matrixes of linear indices (row-wise) and whose second plane shows repeated indices for each
% n·m 'element'. Something like:
%
%   1 2 .. m | 1 2 .. m | ... | 1 2 .. m       1 1 ...  1 | 2 2 ...  2 | ... | M M .. M
%   m+1 .... | m+1 .... | ... | m+1 ....       1 ...      | 2 ...      | ... | M ...
%   ....  nm | ....  nm | ... | ....  nm       ...      1 | ...      2 | ... | ...    M
%   ------------------------------------       ----------------------------------------
%      ...   |    ...   | ... |   ...              ...    |     ...    | ... |   ...
%   ------------------------------------       ----------------------------------------
%   1 2 .. m | 1 2 .. m | ... | 1 2 .. m       (N-1)M+1.. | (N-1)M+2.. | ... | NM .. NM
%   m+1 .... | m+1 .... | ... | m+1 ....       (N-1)M+1.. | (N-1)M+2.. | ... | NM .. NM
%   ....  nm | ....  nm | ... | ....  nm       (N-1)M+1.. | (N-1)M+2.. | ... | NM .. NM

    Nn = N*n; Mm = M*m;

    % Element indices inside each group
    block = repmat(reshape(1:m*n,m,n)',N,M);

    yy = floor((0:Nn-1)/n+1);   % row indices
    xx = floor((0:Mm-1)/m+1);   % column indices
    [xx,yy] = meshgrid(xx,yy);

    % Group indices
    block(:,:,2) = reshape(sub2ind([M,N],xx(:),yy(:)),[Nn,Mm]);
end
