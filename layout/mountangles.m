function [surftilt,surfaz] = mountangles(Trck,sunaz,sunel)
% [SURFTILT,SURFAZ] = MOUNTANGLES(TRCK,SUNAZ,SUNEL) - Wrapper for MOUNTROTATIONS, returns
%   [Nt*,Nm*] arrays of surface-normal tilt(s) and azimuth(s).

    R = mountrotations(Trck,sunaz,sunel);
    surftilt = permute(atan2d(hypot(R(1,3,:,:),R(2,3,:,:)),R(3,3,:,:)),[4,3,1,2]);
    surfaz = permute(atan2d(R(2,3,:,:),R(1,3,:,:)),[4,3,1,2]);
end