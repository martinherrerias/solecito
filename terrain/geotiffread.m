function [I]=geotiffread(filename)
% GEOTIFF_READ: read geotiff using imread and assign map info from infinfo.

% output:
% I.z, image data
% I.x, x coordinate in map
% I.y, y coordinate in map
% I.info, misc. info

% imshow(I.z, 'xdata', I.x, 'ydata', I.y);
% shows image with map coordinate

% Version by Yushin Ahn, ahn.74@osu.edu
% Glacier Dynamics Laboratory, 
% Byrd Polar Resear Center, Ohio State University 
% Referenced enviread.m (Ian Howat)

    Tinfo = imfinfo(filename);
    info.samples = Tinfo(1).Width;
    info.lines = Tinfo(1).Height;
    info.imsize = Tinfo(1).Offset;
    info.bands = Tinfo(1).SamplesPerPixel;

    sub = [1, info.samples, 1, info.lines];

    info.map_info.dx = Tinfo(1).ModelPixelScaleTag(1);
    info.map_info.dy = Tinfo(1).ModelPixelScaleTag(2);
    info.map_info.mapx = Tinfo(1).ModelTiepointTag(4);
    info.map_info.mapy = Tinfo(1).ModelTiepointTag(5);

    xm = info.map_info.mapx;
    ym = info.map_info.mapy;
    xx = xm + ((0:info.samples-1).*info.map_info.dx);
    yy = ym - ((0:info.lines  -1).*info.map_info.dy);

    I.x = xx(sub(1):sub(2));
    I.y = yy(sub(3):sub(4));
    I.z = imread(filename);
    I.info = info;
end





