function D = randomdefects(NmNp,NeNc,P,AB)
% D = RANDOMDEFECTS(NmNp,NeNc,P)
% Generate a random boolean NmNp by NeNc sparse matrix D representing the prevalence of 
% cell-defects in NmNp modules with NeNc cells each, with probability (by module) P of
% having one or more defective cells. That is:
%   D ~ Bernoulli(1-(1-P)^(1/NeNc)) such that nnz(any(D,2))/NmNp ~ P
%
% NmNp and NeNc can be positive integers, or 2-vectors, i.e. [Nm,Np], [Ne,Nc].
%
% D = RANDOMDEFECTS(NmNp,NeNc,P,AB)
% Return instead a sparse-double matrix, where non-zero elements of the boolean D are 
% replaced by Beta(A,B)-distributed random 'defect-strength' values. 
% AB must be a 2-vector [A,B] with the distribution coefficients. e.g. [1,1] for uniform dist.
%
% EXAMPLE: given 500 mounts with 20 modules, 72-cells each, simulate ~5% of modules with Beta(2,5) 
% distributed cell-defects:
%   D = randomdefects([500 20],[3,24],0.05,[2,5]);   % 10000x72 sparse double matrix
%   nnz(any(D,2))/size(D,1)                          % actual percentage of defective modules
%   hist(nonzeros(D));                               % distribution of defect strength

    narginchk(3,4);
    assert(isscalar(P) && P > 0 && P < 1,'Bad P: scalar probability required');
    oksize = @(x,n) all(mod(x,1)==0) && all(x > 0) && isscalar(x) || (isvector(x) && numel(x) == 2);
    assert(oksize(NmNp),'Bad NmNp: expecting integer positive scalar or 2-vector [Nm,Np]');
    assert(oksize(NeNc),'Bad NeNc: expecting integer positive scalar or 2-vector [Ne,Nc]');
    
    if nargin > 3 && ~isempty(AB)
        assert(isvector(AB) && all(AB > 0),'Bad AB: non-negative beta-params required');
    else
        AB = [];
    end
    
    sz = [prod(NmNp),prod(NeNc)];
    p = 1-(1-P)^(1/prod(NeNc));
    D = sparse(rand(sz) < p);
    
    if ~isempty(AB)
        idx = find(D);
        D = D*1; % sparse double
        D(idx) = betarnd(AB(1),AB(2),[numel(idx),1]);
    end
end