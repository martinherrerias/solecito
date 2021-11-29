function [R,C,Xb] = MinBoundCircle(X)
% [R,C,XB] = MINBOUNDCIRCLE(X) - Compute exact minimum bounding circle of a 2D point cloud using 
%   Welzl's algorithm [1]. 
%
% INPUT:
%   - X     : M-by-2 list of point co-ordinates, where M is the total number of points.
%
% OUTPUT:
%   - R     : radius of the minimum bounding circle of X.
%   - C     : 1-by-2 vector specifying centroid co-ordinates of the minimum bounding circle of X.
%   - Xb    : subset of X, listing K-by-2 list of point co-ordinates from which R and C were 
%             computed.
%
% REREFERENCES:
% [1] Welzl, E. (1991), 'Smallest enclosing disks (balls and ellipsoids)',
%     Lecture Notes in Computer Science, Vol. 555, pp. 359-370
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% MIT License, Copyright (c) 2019 Anton Semechko

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.


    if nargin<1 || isempty(X) || ~isnumeric(X) || ~ismatrix(X) 
        error('Unrecognized format for 1st input argument (X)')
    elseif size(X,2)~=2
        error('This function only works for 2D data')
    elseif sum(~isfinite(X(:)))>0
        error('Point data contains NaN and/or Inf entries. Remove them and try again.')    
    end

    % Get minimum bounding circle
    % -------------------------------------------------------------------------

    if size(X,1)>2

        % Check points for collinearity
        Xo=mean(X,1);
        dX=bsxfun(@minus,X,Xo);
        [U,D]=svd((dX'*dX)/size(X,1),0);
        D=diag(D);
        if D(2)<1E-15
            dx=dX*U(:,1);
            [dx_min,id_min]=min(dx);
            [dx_max,id_max]=max(dx);
            R=(dx_max-dx_min)/2;
            C=U(:,1)'*(dx_max+dx_min)/2 + Xo;
            Xb=X([id_min;id_max],:);
            return
        end

        % Get convex hull of the point set
        F=convhull(X);
        F=unique(F(:));
        X=X(F,:);

    end

    try
        % Remove duplicates
        X=uniquetol(X,eps,'ByRows',true);
    catch
        % older version of Matlab; 'uniquetol' is unavailable
    end

    if size(X,1)<3
        [R,C]=FitCircle2Points(X); 
        Xb=X;   
        return
    end

    % Randomly permute the point set
    idx=randperm(size(X,1));
    X=X(idx(:),:);
    if size(X,1)<1E3
        try

            % Center and radius of the circle
            [R,C]=B_MinCircle(X,[]);

            % Co-ordinates of the points used to compute parameters of the 
            % minimum bounding circle
            D=sum(bsxfun(@minus,X,C).^2,2);
            [D,idx]=sort(abs(D-R^2));
            Xb=X(idx(1:4),:);
            D=D(1:4);
            Xb=Xb(D<1E-6,:);
            [~,idx]=sort(Xb(:,1));
            Xb=Xb(idx,:);
            return
        catch
        end
    end

    % If we got to this point, then either size(X,1)>=1E3 or recursion depth 
    % limit was reached. So need to break-up point-set into smaller sets and 
    % then recombine the results.
    M=size(X,1);
    dM=max(min(floor(M/4),300),3);
    res=mod(M,dM);
    n=ceil(M/dM);  
    idx=dM*ones(1,n);
    if res>0
        idx(end)=res;
    end

    if res<0.25*dM && res>0
        idx(n-1)=idx(n-1)+idx(n);
        idx(n)=[];
        n=n-1;
    end

    X=mat2cell(X,idx,2);
    Xb=[];
    for i=1:n

        % Center and radius of the circle
        [R,C,Xi]=B_MinCircle([Xb;X{i}],[]);    

        % 40 points closest to the circle
        if i<1
            D=abs(sum(bsxfun(@minus,Xi,C).^2,2)-R^2);
        else
            D=abs(sqrt(sum(bsxfun(@minus,Xi,C).^2,2))-R);
        end
        [D,idx]=sort(D);
        Xb=Xi(idx(1:min(40,numel(D))),:);

    end
    D=D(1:3);
    Xb=Xb(D/R*100<1E-3,:);
    [~,idx]=sort(Xb(:,1));
    Xb=Xb(idx,:);
end

function [R,C,P]=B_MinCircle(P,B)

    if size(B,1)==3 || isempty(P)
        [R,C]=FitCircle2Points(B); % fit circle to boundary points
        return
    end

    % Remove the last (i.e., end) point, p, from the list
    P_new=P;
    P_new(end,:)=[];
    p=P(end,:);

    % Check if p is on or inside the bounding circle. If not, it must be
    % part of the new boundary.
    [R,C,P_new]=B_MinCircle(P_new,B); 
    if isnan(R) || isinf(R) || R<=eps
        chk=true;
    else
        chk=norm(p-C)>(R+eps);
    end

    if chk
        B=[p;B];
        [R,C]=B_MinCircle(P_new,B);
        P=[p;P_new];
    end
end

function [R,C] = FitCircle2Points(X)
% Fit a circle to a set of 2 or at most 3 points in 3D space. Note that
% point configurations with 3 collinear points do not have well-defined 
% solutions (i.e., they lie on circles with infinite radius).
%
% INPUT:
%   - X     : M-by-2 array of point coordinates, where M<=3.
%
% OUTPUT:
%   - R     : radius of the circle. R=Inf when the circle is undefined, as 
%             specified above.
%   - C     : coordinates of the circle centroid. C=nan(1,2) when the 
%             circle is undefined, as specified above.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)


    N=size(X,1);
    if N>3
        error('Input must a N-by-2 array of point coordinates, with N<=3')
    end

    % Empty set
    if isempty(X)
        C=nan(1,2);
        R=nan; 
        return
    end

    % A single point
    if N==1
        C=X;
        R=0;
        return
    end

    % Line segment
    if N==2
        C=mean(X,1);
        R=norm(X(2,:)-X(1,:))/2;
        return
    end

    % Remove duplicate vertices, if there are any
    D=bsxfun(@minus,permute(X,[1 3 2]),permute(X,[3 1 2]));
    D=sqrt(sum(D.^2,3));
    D(1:(N+1):end)=Inf;
    chk=D<=1E-12;
    if sum(chk(:))>0
        for i=1:(N-1)
            if size(X,1)<=i, break; end
            idx=chk(i,:);
            idx(1:i)=false;
            idx=find(idx);
            chk(idx,:)=[];
            chk(:,idx)=[];
            X(idx,:)=[];
        end
        [R,C]=FitCircle2Points(X);
        return
    end


    % Three unique, though possibly collinear points
    tol=1E-2; % collinearity threshold (in degrees)

    % Check for collinearity
    D12=X(2,:)-X(1,:); D12=D12/norm(D12);
    D13=X(3,:)-X(1,:); D13=D13/norm(D13);

    chk=abs(D12*D13(:));
    chk(chk>1)=1;
    if acos(chk)*(180/pi)<tol
        R=inf;
        C=nan(1,2);
        return
    end

    % Circle centroid
    A=2*bsxfun(@minus,X(2:3,:),X(1,:));
    b=sum(bsxfun(@minus,X(2:3,:).^2,X(1,:).^2),2);
    C=(A\b)';

    % Circle radius
    R=sqrt(sum((X(1,:)-C).^2,2));
end
