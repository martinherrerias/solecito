
N = 50;   % points per curve
M = 1000;   % number of curves

X = sort(rand(N,M),1);
Y = sort(rand(N,M),1,'descend');
C = cell(M,1);
P = mdpwl(); P(M,1) = mdpwl();


fprintf('\n\n\n\n');
fprintf('griddedInterpolant(),..: ');
tic()
for j = 1:M
    C{j} = griddedInterpolant();
    C{j}.GridVectors = {X(:,j)};
    C{j}.Values = Y(:,j);
end
toc()
fprintf('             evaluation: ');
tic()
for j = 1:M
    z = C{j}(rand(N,1));
end
toc()


fprintf('\n');
fprintf('griddedInterpolant(x,y): ');
tic()
for j = 1:M
    C{j} = griddedInterpolant(X(:,j),Y(:,j));
end
toc()
fprintf('             evaluation: ');
tic()
for j = 1:M
    z = C{j}(rand(N,1));
end
toc()


fprintf('\n');
fprintf('    interp1(x,y,pp),...: ');
tic()
for j = 1:M
    C{j} = interp1(X(:,j),Y(:,j),'linear','pp'); %#ok
end
toc()
fprintf('             evaluation: ');
tic()
for j = 1:M
    ppval(C{j},rand(N,1));
end
toc()

path = strrep(fileparts(which('mdpwl')),'simulation/07','simulation (backup-unused)/07');
cleaner = onCleanup(@() rmpath(path)); addpath(path);

fprintf('\n');
fprintf('            mkivpp(x,y): ');
tic()
for j = 1:M
    C{j} = mkivpp(X(:,j),Y(:,j),0);
end
toc()
fprintf('             ppval(...): ');
tic()
for j = 1:M
    ppval(C{j},rand(N,1));
end
toc()

fprintf('\n');
fprintf('             mdpwl(x,y): ');
tic()
for j = 1:M
    P(j) = mdpwl(X(:,j),Y(:,j),0);
end
toc()
fprintf('         mdpwl.val(...): ');
tic()
for j = 1:M
    P(j).val(rand(N,1));
end
toc()

% fprintf('\n');
% fprintf('         mex_mdpwl(x,y): ');
% tic()
% for j = 1:M
%     P(j) = mex_mdpwl(X(:,j),Y(:,j),0);
% end
% toc()
% fprintf('         mdpwl.val(...): ');
% tic()
% for j = 1:M
%     P(j).val(rand(N,1));
% end
% toc()



