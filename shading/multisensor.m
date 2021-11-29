function [n,P] = multisensor(varargin)
% [N,P] = MULTISENSOR(Za,Aa,Zb,Ab,Zc,Ac...) - return normal vectors N and Voronoi polygons P
%   built around those normals, for a set of points with zenith angle(s) [Za,Zb,Zc,..] and azimuth
%   angle(s) [Aa,Ab,Ac,..].
%
% INPUT: 
%   Za,Aa,Zb,Ab ... - a series of PAIRED zenith/azimuth angles (degrees) for any number of sensors.
%       Each Zx,Ax can be either a scalar (that will be expanded to match its pair) or a vector.
%
% OUTPUT:
%   n - 3xm array of sensor normal vectors (pointing outwards).
%   P - mx1 cell array of polygon3d objects (spherical Voronoi polygons)
%
% EXAMPLE: [n,P] = MULTISENSOR(0,0,90,0:90:270) - returns the 5 cardinal points (including zenith).

    okargs = @(x) isvector(x) && isnumeric(x) && all(isfinite(x));
    assert(mod(nargin,2) == 0 && all(cellfun(okargs,varargin)),...
        'Expecting arguments as pairs of angle vectors');
    for j = 1:2:nargin
        if isscalar(varargin{j})
            varargin{j}(1:numel(varargin{j+1})) =  varargin{j};
        elseif isscalar(varargin{j+1})
            varargin{j+1}(1:numel(varargin{j})) =  varargin{j+1};
        end
    end
    varargin = cellfun(@(x) reshape(x,1,[]),varargin,'unif',0);
    z = cat(2,varargin{1:2:end});
    az = cat(2,varargin{2:2:end});

    [x,y,z] = sph2cart(az*pi/180,pi/2-z*pi/180,1);
    
    n = [x;y;z];
    [~,~,P] = voronoisphere(n,'resolution',1.5); % 1.5 rad just to split full-wedges
    P = cellfun(@polygon3d,P','unif',0);
    %P = mat2cell(P(:),ones(numel(P),1));
end