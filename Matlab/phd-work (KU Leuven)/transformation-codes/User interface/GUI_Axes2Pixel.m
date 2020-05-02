function p = GUI_Axes2Pixel(ax, v, vector)
%GUI_Axes2Pixel  Transform vertices or vectors from axes to pixel space.
%
%   Syntax:
%    p = GUI_Axes2Pixel(ax, v, vector)
%
%   Input:
%    ax:     Axes handle.
%    v:      N-by-3 array containing vector or vertex coordinates in ax's
%            data space. Each row represents a point; the first, second
%            and third columns represent X-, Y- and Z-coordinates
%            respectively.
%    vector: Logical indicating whether the elements of v are vectors. If
%            set to true, the coordinate origin will not be translated;
%            otherwise, the data coordinates will be treated as vertices.
%            Optional, defaults to false.
%
%   Output:
%    p: N-by-3 array containing pixel coordinates. Each row represents a
%       point; the first, second and third columns represent X-, Y- and
%       Z-coordinates respectively. On the screen, X corresponds with
%       right, Y corresponds with up and Z points out of the screen towards
%       the viewer.
%
%   Effect: This function will transform the provided data coordinates in
%   axes ax into data coordinates. Points on the camera view plane will
%   transform into pixels with Z-coordinates equal to 0. This function is
%   the opposite of GUI_Pixel2Axes.
%
%   Example:
%    GUI_Axes2Pixel(gca, get(gca, 'CurrentPoint'))
%
%   Dependencies: GUI_Pixel2AxesTransform.m
%
%   Known parents: none

%Created on 23/05/2006 by Ward Bartels.
%WB, 16/07/2007: Rewrite based on GUI_Pixel2AxesTransform.m.
%Stabile, fully functional.


%Set default for vector
vector = nargin>=3 && vector;

%Obtain rotation matrix and apply translation <<GUI_Pixel2AxesTransform.m>>
if vector
    R = GUI_Pixel2AxesTransform(ax);
else
    [R, T] = GUI_Pixel2AxesTransform(ax);
    v = v-T(ones(size(v, 1), 1),:);
end

%Apply rotation
p = v/R;


% %Original code (wb20070716), begin:
% 
% %Get needed transforms (this is not supported by TMW)
% xform = get(ax,'x_RenderTransform'); 
% offset = get(ax,'x_RenderOffset');
% scale = get(ax,'x_RenderScale');
% 
% %Equivalent: nvert = vert/scale-offset;
% nvert(:,1) = V(:,1)./scale(1)-offset(1);
% nvert(:,2) = V(:,2)./scale(2)-offset(2);
% nvert(:,3) = V(:,3)./scale(3)-offset(3);
% 
% %Equivalent: xvert = xform*xvert;
% w = xform(4,1)*nvert(:,1)+xform(4,2)*nvert(:,2)+xform(4,3)*nvert(:,3)+xform(4,4);
% xvert(:,1) = xform(1,1)*nvert(:,1)+xform(1,2)*nvert(:,2)+xform(1,3)*nvert(:,3)+xform(1,4);
% xvert(:,2) = xform(2,1)*nvert(:,1)+xform(2,2)*nvert(:,2)+xform(2,3)*nvert(:,3)+xform(2,4);
% 
% %Fail-safe (w may be 0 for perspective plots )
% ind = w==0;
% w(ind) = 1;
% xvert(ind,:) = 0;
% 
% %Eliminate scaling coordinate
% p(:,1) = xvert(:,1)./w;
% p(:,2) = xvert(:,2)./w;
% 
% %Original code (wb20070713), end.