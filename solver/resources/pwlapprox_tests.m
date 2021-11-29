
%[y,dy] = circle(1-eps)
p = pwlapprox(@circle,0,1-eps,1e-6);
plot(p,'r.')

% function [y,m] = circle(x)
%     y = sqrt(1-x.^2);
%     m = -x./y;
% end

function y = circle(x)
    y = sqrt(1-x.^2);
    %m = -x./y;
end