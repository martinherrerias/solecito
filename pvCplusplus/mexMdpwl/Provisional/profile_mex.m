N = 100;
M = 1;

X1 = rand(N,M);
X1 = sort(X1,1);
Y1 = rand(N,M);
Y1 = sort(Y1,1,'descend');

X2 = rand(N, M);
X2 = sort(X2,1);
Y2 = rand(N,M);
Y2 = sort(Y2,1,'descend');

X3 = rand(N, M);
X3 = sort(X3,1);
Y3 = rand(N,M);
Y3 = sort(Y2,1,'descend');

X4 = rand(N, M);
X4 = sort(X4,1);
Y4 = rand(N,M);
Y4 = sort(Y4,1,'descend');

X5 = rand(N, M);
X5 = sort(X5,1);
Y5 = rand(N,M);
Y5 = sort(Y5,1,'descend');

X6 = rand(N, M);
X6 = sort(X6,1);
Y6 = rand(N,M);
Y6 = sort(Y6,1,'descend');

X7 = rand(N, M);
X7 = sort(X7,1);
Y7 = rand(N,M);
Y7 = sort(Y7,1,'descend');


X8 = rand(N, M);
X8 = sort(X8,1);
Y8 = rand(N,M);
Y8 = sort(Y8,1,'descend');


A = cell(M,1);
B = cell(M,1);
C = cell(M,1);
D = cell(M,1);
E = cell(M,1);
F = cell(M,1);
G = cell(M,1);
H = cell(M,1);


% fprintf('\ngriddedInterpolant(x,y): ');
% tic()
% for j = 1:M
%     B{j} = griddedInterpolant(X(:,j),Y(:,j));
% end
% toc()
% 
% fprintf('griddedInterpolant(),...: ');
% tic()
% for j = 1:M
%     A{j} = griddedInterpolant();
%     A{j}.GridVectors = {X(:,j)};
%     A{j}.Values = Y(:,j);
% end
% toc()
% 
% fprintf('interp1(x,y,pp),...: ');
% tic()
% for j = 1:M
%     C{j} = interp1(X(:,j),Y(:,j),'linear','pp'); %#ok
% end
% toc()
% 
% fprintf('eval: ');
% tic()
% for j = 1:M
%     z = A{j}(rand(N,1));
% end
% toc()

%fprintf('\nmkivpp(x,y):');
%tic()
%for j = 1:M
%    A{j} = mkivpp(X(:,j),Y(:,j));
%end
%toc()


%fprintf('mdpwl(x,y):');

%for j = 1:M
%    A{j} = mexaddparallel([X1(:,j),Y1(:,j)]);
%    B{j} = mexaddseries([X1(:,j),Y1(:,j)], -4, 4);
%end

tic()
for j = 1:M
    
    [out1, out2, size] = mexaddparallel(8, [X1(:,j),Y1(:,j)], [X2(:,j),Y2(:,j)], [X3(:,j),Y3(:,j)],[X4(:,j),Y4(:,j)],[X5(:,j),Y5(:,j)],[X6(:,j),Y6(:,j)],[X7(:,j),Y7(:,j)],[X8(:,j),Y8(:,j)], [0, 0, 0], [-0.6, 0.7, -Inf, Inf]);
%    Coutput = mdpwl(out1(1:size)', out2(1:size)', 0);
end
toc()

disp(size);

tic()
for j = 1:M
	A{j} = mdpwl(X1(:,j),Y1(:,j));
    B{j} = mdpwl(X2(:,j),Y2(:,j));
    C{j} = mdpwl(X3(:,j),Y3(:,j));
    D{j} = mdpwl(X4(:,j),Y4(:,j)); 
    E{j} = mdpwl(X5(:,j),Y5(:,j));
    F{j} = mdpwl(X6(:,j),Y6(:,j));
    G{j} = mdpwl(X7(:,j),Y7(:,j));
    H{j} = mdpwl(X8(:,j),Y8(:,j));

    matlaboutput = addparallel([A{j},B{j},C{j},D{j},E{j},F{j},G{j},H{j}], [-0.6, 0.7, -Inf, Inf], [0, 0, 0]);
end
toc()
disp(matlaboutput.n);

%fprintf('mex(x,y):');
%tic()
%for j = 1:M
%    B{j} = mexdouglaspeucker([X(:,j),Y(:,j)]);
%    B{j} = mdpwl(X(:,j),Y(:,j),0);
%end
%toc()

% fprintf('mex_mdpwl(x,y):');
% tic()
% for j = 1:M
%     B{j} = mex_mdpwl(X(:,j),Y(:,j),0);
% end
% toc()

%fprintf('\nppval():');
%tic()
%for j = 1:M
%    ppval(A{j},rand(N,1));
%end
%toc()

%fprintf('mdpwl.val():');
%tic()
%for j = 1:M
%    B{j}.val(rand(N,1));
%end
%toc()


