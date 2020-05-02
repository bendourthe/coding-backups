function center = GUI_CameraCenterPixel(h)
%GUI_CameraCenterPixel  Camera center point in figure pixel coordinates.
%wb20070716
%
%   Syntax:
%    center = GUI_CameraCenterPixel(ax)
%
%   Input:
%    ax: Axes handle.
%
%   Output:
%    center: 2-element row vector containing the camera center point on the
%            figure in pixel coordinates.
%
%   Effect: This function will calculate the camera center point in pixel
%   coordinates, in the figure's coordinate system. Both the camera
%   position and the camera target are projected onto this point.
%
%   Dependencies: none
%
%   Known parents: GUI_Pixel2AxesTransform.m
%                  GUI_FreeRotate.m

%Created on 16/07/2007 by Ward Bartels.
%Stabile, fully functional.


%Get figure handle
fig = ancestor(h, 'figure');

%Get positions of axes and its parents
pos = zeros(0, 4);
while h~=fig
    h_par = get(h, 'Parent');
    pos = [pos; hgconvertunits(fig, get(h, 'Position'), get(h, 'Units'), 'pixels', h_par)];
    h = h_par;
    % units = get(h, 'Units');
    % set(h, 'Units', 'pixels');
    % pos = [pos; get(h, 'Position')];
    % set(h, 'Units', units);
    % h = get(h, 'Parent');
end

%Calculate camera center point in figure pixel coordinates
center = sum(pos(:,[1 2]), 1)+pos(1,[3 4])/2-[0 1];