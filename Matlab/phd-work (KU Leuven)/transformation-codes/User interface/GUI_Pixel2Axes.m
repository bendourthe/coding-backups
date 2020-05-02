function v = GUI_Pixel2Axes(ax, p, vector)
%GUI_Pixel2Axes  Transform vertices or vectors from pixel to axes space.
%
%   Syntax:
%    v = GUI_Pixel2Axes(ax, p, vector)
%
%   Input:
%    ax:     Axes handle.
%    p:      N-by-3 array containing pixel coordinates. Each row represents a
%            point; the first, second and third columns represent X-, Y- and
%            Z-coordinates respectively. On the screen, X corresponds
%            with right, Y corresponds with up and Z points out of the
%            screen towards the viewer. The last column may be omitted,
%            which will result in Z being set to 0 for all points.
%    vector: Logical indicating whether the elements of p are vectors. If
%            set to true, the coordinate origin will not be translated;
%            otherwise, the pixel coordinates will be treated as vertices.
%            Optional, defaults to false.
%
%   Output:
%    v: N-by-3 array containing vector or vertex coordinates in ax's data
%       space. Each row represents a point; the first, second and third
%       columns represent X-, Y- and Z-coordinates respectively.
%
%   Effect: This function will transform the provided pixel coordinates
%   into data coordinates in axes ax. Pixels with Z-coordinates equal to 0
%   will be transformed onto the camera view plane. This function is the
%   opposite of GUI_Axes2Pixel.
%
%   Example:
%    GUI_Pixel2Axes(gca, get(gcf, 'CurrentPoint'))
%
%   Dependencies: GUI_Pixel2AxesTransform.m
%
%   Known parents: GUI_ManipulateObjects.m

%Created on 16/07/2007 by Ward Bartels.
%Stabile, fully functional.


%Set default for vector
vector = nargin>=3 && vector;

%Obtain transformation matrices <<GUI_Pixel2AxesTransform.m>>
if vector
    R = GUI_Pixel2AxesTransform(ax);
else
    [R, T] = GUI_Pixel2AxesTransform(ax);
end

%Apply rotation
v = p*R(1:size(p, 2),:);

%Apply translation if necessary
if ~vector
    v = v+T(ones(size(v, 1), 1),:);
end