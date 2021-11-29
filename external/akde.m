function [pdf,X1,X2] = akde(X,grid,gam)
% [pdf,X1,X2] = AKDE(X,grid,gam) - adaptive kernel density estimation in high dimensions;
%   optimal accuracy/speed tradeoff, controlled via parameter "gam"
%
%   Taken, with minimum modifications (readability, memory footprint) from [1].
%
% INPUTS:   X  - n x d matrix
%
%         grid - 'm' points of dimension 'd' over which pdf is computed;
%                default provided only for 2-dimensional data;
%                see example on how to construct it in higher dimensions
%
%          gam - cost/accuracy tradeoff parameter, where gam < n;
%                default value is gam = ceil(n^(1/2)); larger values
%                may result in better accuracy, but always reduce speed;
%                to speedup the code, reduce the value of "gam"; 
%
% OUTPUT: pdf   - the value of the estimated density at 'grid'
%         X1,X2 - grid only for 2 dimensional data
%
% EXAMPLES:
%   % in 2 dimensions:
%       L = chol([1,-0.999;-0.999,1],'lower');
%       L1=chol([1,0.999;0.999,1],'lower');
%       data = [(L1*randn(10^3,2)')';
%       (L*randn(10^3,2)')'*2;rand(10^4,2)*5-2.5];
%       [pdf,X1,X2] = akde(data);
%       pdf=reshape(pdf,size(X1));
%       contour(X1,X2,pdf,20)
%
%  % in 3 dimensions:
%      data = [randn(10^3,3);
%      randn(10^3,3)/2+2]; % three dimensional data
%      [~,d] = size(data); 
%      ng = 100; % total grid points = ng^d
%      MAX = max(data,[],1); 
%      MIN = min(data,[],1); 
%      scaling = MAX-MIN;
%      % create meshgrid in 3-dimensions
%      [X1,X2,X3] = meshgrid(MIN(1):scaling(1)/(ng-1):MAX(1),...
%          MIN(2):scaling(2)/(ng-1):MAX(2),MIN(3):scaling(3)/(ng-1):MAX(3));
% 
%      grid = reshape([X1(:),X2(:),X3(:)],ng^d,d); % create points for plotting
%      pdf = akde(data,grid); % run adaptive kde
%      pdf = reshape(pdf,size(X1)); % reshape pdf for use with meshgrid
%      for iso = 0.005:0.005:0.015 
%      % isosurfaces with pdf = 0.005,0.01,0.015
%          isosurface(X1,X2,X3,pdf,iso),view(3),alpha(.3),box on,hold on
%          colormap cool
%      end
%
% REF:
%  [1] Zdravko Botev (2021). Kernel Density Estimator for High Dimensions 
%      https://www.mathworks.com/matlabcentral/fileexchange/58312-kernel-density-estimator-
%      for-high-dimensions, MATLAB Central File Exchange. Retrieved October 2, 2021. 
%  [2] Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010), "Kernel density estimation via 
%      diffusion, Annals of Statistics, Volume 38, Number 5, pages 2916-2957.
%
% Copyright (c) 2016, Zdravko Botev
% Copyright (c) 2015, Zdravko Botev
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% * Neither the name of University of New South Wales nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% * Neither the name of The University of New South Wales nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    [n,d] = size(X);
    
    % begin scaling preprocessing
    MAX = max(X,[],1);
    MIN = min(X,[],1);
    scaling = MAX-MIN;
    MAX = MAX+scaling/10;
    MIN = MIN-scaling/10;
    scaling = MAX-MIN;
    X = X-MIN;
    X = X./scaling;
    
    if (nargin<2) || isempty(grid) % failing to provide grid
        warning('Assuming data is 2 dimensional. For higher dimensions, provide a grid as in example.')
        % create meshgrid in 2-dimensions
        [X1,X2]=meshgrid(MIN(1):scaling(1)/(2^7-1):MAX(1),...
               MIN(2):scaling(2)/(2^7-1):MAX(2));
        grid=reshape([X1(:),X2(:)],2^14,d); % create grid for plotting
    end
    grid = (grid - MIN)./scaling;
    if nargin<3 % failing to provide speed/accuracy tradeoff
        gam=ceil(n^(1/2));
    end
    
    % end preprocessing
    % algorithm initialization
    del = 0.1/n^(d/(d+4));
    perm = randperm(n);
    mu = X(perm(1:gam),:);
    w = rand(1,gam);
    w = w/sum(w);
    Sig = rand(d,d,gam).*eye(d)*del;
    ent=-Inf;
    
    fprintf('----------------------------\n');
    fprintf('Iter.    Tol.      Bandwidth \n');
    for iter=1:1500 % begin algorithm
        Eold=ent;
        [w,mu,Sig,del,ent]=regEM(w,mu,Sig,del,X); % update parameters
        err=abs((ent-Eold)/ent); % stopping condition
        fprintf('%4i    %8.2e   %8.2e\n',iter,err,del);
        if (err<10^-4) || (iter>200), break, end
    end
    
    fprintf('----------------------------\n');
    % now output density values at grid
    pdf = probfun(grid,w,mu,Sig)/prod(scaling); % evaluate density
    % del=del*scaling; % adjust bandwidth for scaling
end

function pdf=probfun(x,w,mu,Sig)
    [gam,d] = size(mu);
    pdf = 0;
    for k=1:gam
        L = chol(Sig(:,:,k));
        s = diag(L);
        logpdf = -0.5*sum(( (x - mu(k,:))/L ).^2,2) + log(w(k)) - sum(log(s)) - d*log(2*pi)/2;
        pdf = pdf + exp(logpdf);
    end
end

function [w,mu,Sig,del,ent]=regEM(w,mu,Sig,del,X)
    [gam,~] = size(mu);
    [n,d] = size(X);

    maxll = nan(n,1);
    maxlsig = nan(n,1);
    psigd = zeros(n,1);
    p = zeros(n,gam);

    for i = 1:gam
        L = chol(Sig(:,:,i));
        Xcentered = X - mu(i,:);
        xRinv = Xcentered /L; 
        xSig = sum((xRinv /L').^2,2)+eps;

        log_lh_i = -0.5*sum(xRinv.^2, 2)-sum(log(diag(L)))...
            +log(w(i))-d*log(2*pi)/2-.5*del^2*trace((eye(d)/L)/L');
        log_sig_i = log_lh_i+log(xSig);

        maxll = max(log_lh_i,maxll);
        maxlsig = max(log_sig_i,maxlsig);

        p(:,i) = log_lh_i;
        psigd = psigd + exp(log_sig_i);
    end
    p = exp(p - maxll); 

    psigd = exp(log(psigd)-maxlsig);

    density = sum(p,2);  
    logpdf = log(density) + maxll;
    logpsigd = log(psigd) + maxlsig;
    p = p./density;
    ent = sum(logpdf);
    w = sum(p,1);

    for i = find(w>0)
        mu(i,:) = p(:,i)'*X/w(i);  %compute mu's
        Xcentered = X - mu(i,:);
        Xcentered = sqrt(p(:,i)).*Xcentered;
        Sig(:,:,i)=Xcentered'*Xcentered/w(i)+del^2*eye(d); % compute sigmas;
    end
    w = w/sum(w);
    curv = mean(exp(logpsigd-logpdf)); % estimate curvature
    del = 1/(4*n*(4*pi)^(d/2)*curv)^(1/(d+2));
end