function [V,T] = synterrain(d0, h0, r0, rr, n0, m, n, Kd, Kh, cdf)
% [x, y, h] = SYNTERRAIN(d0, h0, r0, rr, n0, m, n, Kd, Kh)
% Generates a TIN (triangulated irragular network) of points that resembles terrain, using a 
% variation of fractional-Brownian-Noise.
% The function is based on the GENERATE_TERRAIN algoritm of Tucker McClure (2012) [1] but offers 
% a full range of parameters for higher output control.
%
% Inputs:
% n         - Number of iterations of algorithm to perform, akin to the
%             order of magnitudes represented by the variation. Anything
%             beyond 7 will produce detail too fine to notice, while the
%             algorithm will require ~3 times longer for every iteration.
% mesh_size - Size of the output mesh (e.g., 512 for 512-by-512)
% h0        - Initial elevation
% r0        - Initial roughness (how much terrain can vary in a step)
% rr        - Roughness roughness (how much roughness can vary in a step)
%
% For an example, see terrain_generation_introduction.m or 
% html/terrain_generation_introduction.html.
%
% Outputs:
% x, y, and h    - Vectors of points comprising terrain
% xm, ym, and hm - Meshes over landscape, useful for surf(...)
%
% [1] https://www.mathworks.com/matlabcentral/fileexchange/39559-automatic-terrain-generation

    % Set some defaults if the user didn't provide inputs.
    if nargin < 1 || isempty(d0), d0 = 1; end                   % Approx diameter
    if nargin < 2 || isempty(h0), h0 = 0.1 * d0 * rand(); end   % Elevation
    if nargin < 3 || isempty(r0), r0 = h0 * rand(); end         % Roughness
    if nargin < 4 || isempty(rr), rr = r0 * rand(); end         % Roughness roughness
    if nargin < 5 || isempty(n0), n0 = randi(5); end            % Number of initial points
    if nargin < 6 || isempty(m), m = 3; end                     % Children per iteration
    if nargin < 7 || isempty(n), n = 7; end                     % Iterations
    if nargin < 8 || isempty(Kd), Kd = 0.75; end                % Horizontal var. decrease coef.
    if nargin < 9 || isempty(Kh), Kh = 0.5; end                 % Vertical var. decrease coef.
    if nargin < 10, cdf = []; end

    if isempty(cdf), randf = @randn; 
    else, randf = @(m,n) interp1(cdf(:,2),cdf(:,1),rand(m,n));
    end
    
    nf = n0 * (m+1)^n; % Total number of points
    
    % Create initial x, y, and height coordinates and roughness map.
    x = [rand(n0, 1); zeros(nf-n0, 1)];
    y = [rand(n0, 1); zeros(nf-n0, 1)];
    h = [randn(n0, 1); zeros(nf-n0, 1)];
    r = [rr * randn(n0, 1) + r0; zeros(nf-n0, 1)];
    
    % Create new points from old points n times.
    for k = 1:n

        % Calculate the new variance for the x, y random draws and for the h, r random draws.
        dxy = Kd^k;
        dh  = Kh^k;

        % Number of new points to generate
        n_new = m * n0;

        % Parents for new points
        parents = reshape(repmat(1:n0, m, 1), [n_new, 1]);
        
        % Calculate indices for new and existing points.
        new = (n0+1):(n0+n_new);
        old = 1:n0;

        % Generate new x/y values.
        theta  = 2*pi * rand(n_new, 1);
        radius = dxy * (rand(n_new, 1) + 1);
        x(new) = x(parents) + radius .* cos(theta);
        y(new) = y(parents) + radius .* sin(theta);
        
        % Interpolate to find nominal new r and h values and add noise to roughness and height maps.
        r(new) = interpolate(x(old), y(old), r(old), x(new), y(new)) ...
                 + (dh * rr) .* randf(n_new, 1);
        h(new) = interpolate(x(old), y(old), h(old), x(new), y(new)) ...
                 + (dh/dxy) * radius .* r(new) .* randf(n_new, 1);
        
        % Add up the new points.
        n0 = n_new + n0;
    end

    % Normalize the distribution of the points about the median.
    x = (x - median(x))/std(x)*d0;
    y = (y - median(y))/std(y)*d0;
    
    %h = (h - min(h))/std(h)*h0;
    h = (h - median(h))/std(h)*h0;
    
    % Make valleys flat, remove unnecesary vertices
    h(h < 0) = 0;
    T = delaunay(x,y);
    flat = all(h(T) == 0,2);
    keep = false(numel(x),1);
    for j = 1:numel(x)
       inT = any(T == j,2);
       keep(j) = ~all(flat(inT));
    end
    
    V = [x(keep),y(keep),h(keep)];
    T = delaunay(x(keep),y(keep));
end

% Perform our particular type of interpolation.
function vn = interpolate(x0, y0, v0, xn, yn)

    % Make an interpolator for height or roughness. We'll add safe far-away
    % points so we'll always be interpolating and not extrapolating.
    int = scatteredInterpolant([100*[-1 -1 1 1]'; x0], [100*[-1 1 -1 1]'; y0],[zeros(4,1); v0]);

    % Perform the actual interpolation.
    vn = int(xn, yn);
end

