function [R, T] = GUI_Pixel2AxesTransform(ax)
%GUI_Pixel2AxesTransform  Create transformation from pixel to data space.
%
%   Syntax:
%    [R, T] = GUI_Pixel2AxesTransform(ax)
%
%   Input:
%    ax: Axes handle.
%
%   Output:
%    R: 3-by-3 rotation matrix.
%    T: 3-element translation row vector.
%
%   Effect: This function will generate a rotation matrix R and a
%   translation vector T representing a transformation from pixel
%   coordinates to data coordinates in axes ax. A pixel may be represented
%   by a 3-element row vector called p. The three elements are X-, Y- and
%   Z-coordinates. On the screen, X corresponds with left, Y corresponds
%   with up and Z points out of the screen towards the viewer. Data
%   coordinates may be calculated as p*R+T. If the pixel's Z-coordinate is
%   0, it will be transformed onto the camera view plane.
%
%   Note: Due to bugs in MATLAB (incorrect axes CameraViewAngle), this
%   function will return unexpected results when:
%    - The axes' CameraViewAngleMode, DataAspectRatioMode and
%      PlotBoxAspectRatioMode are set to 'auto';
%    - The axes' CameraViewAngleMode is set to 'auto' and
%      PlotBoxAspectRatio is anything other than [1 1 1]
%
%   Dependencies: VectorNorms.m
%                 GUI_CameraCenterPixel.m
%
%   Known parents: GUI_Pixel2Axes.m
%                  GUI_Axes2Pixel.m

%Created on 16/07/2007 by Ward Bartels.
%Stabile, fully functional.


%Store axes DataAspectRatio and CameraPosition
dar = daspect(ax);
camloc = campos(ax);

%Create rotation matrix
R(3,:) = (camloc-camtarget(ax))./dar;
R(1,:) = cross(camup(ax)./dar, R(3,:)); 
R(2,:) = cross(R(3,:), R(1,:));

%Normalize matrix <<VectorNorms.m>>
norms = VectorNorms(R);
R = R./(norms*[1 1 1]);

%Calculate field of view
fov = 2*norms(3)*tand(camva(ax)/2);

%Get axes size in pixels
fig = ancestor(ax, 'figure');
pos = hgconvertunits(fig, get(ax, 'Position'), get(ax, 'Units'), 'pixels', get(ax, 'Parent'));
% units = get(ax, 'Units');
% set(ax, 'Units', 'pixels');
% pos = get(ax, 'Position');
% set(ax, 'Units', units);
pix = min(pos(3), pos(4));

%Scale rotation matrix by field of view
R = R*diag((dar)*fov/pix);

%If translation wasn't requested, return
if nargout<2, return; end

%Calculate translation if requested <<GUI_CameraCenterPixel.m>>
if nargout>=2
    T = camloc-GUI_CameraCenterPixel(ax)*R([1 2],:);
end