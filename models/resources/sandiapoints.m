function [vv,ii] = sandiapoints(pars)
% [VV,II] = SANDIAPOINTS(PARS) - returns the five SANDIA points for the ODM defined by PARS:
% 
%   VV = [ 0    Vmp/2  Vmp   (Vmp+Voc)/2  Voc] 
%   II = [Isc I(Vmp/2) Imp I((Vmp+Voc)/2)  0 ]
%
% PARS is (an array of) structure(s) of adjusted ODM parameters as returned from TRANSLATEODM.
%
% If PARS is an Nx1 vector of structures, VV,II will be Nx5 arrays.
% If PARS is an 1xN vector of structures, VV,II will be 5xN arrays.
% For other, arbitrary sized arrays, VV,II will be [size(PARS),5] arrays.
%
% See also: TRANSLATEODM, ONEDIODEMPP


    vv = zeros(numel(pars),5);
    ii = zeros(numel(pars),5);
    [~,vv(:,3),ii(:,3),vv(:,5),ii(:,1)] = onediodeMPP(pars(:));
    vv(:,2)=vv(:,5)/2;
    vv(:,4)=(vv(:,5)+vv(:,3))/2;
    for j = 1:numel(pars)
        ii(j,[2,4]) = onediodemodel(pars(j),vv(j,[2,4]));
    end
    if j > 1
        if isvector(pars)
            if size(pars,1) == 1, vv = vv'; ii = ii'; 
            else % if size(ODMpar,2) == 1, do nothing
            end
        else
            vv = reshape(vv,[size(pars),5]);
            ii = reshape(ii,[size(pars),5]);
        end
    end
end
