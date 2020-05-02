function [wbdfun, kpfun] = GUI_FreeRotate(ax, enable, coordsys, box)
%GUI_FreeRotate  Enable interactive rotation of a 3-D plot.
%
%   Syntax:
%    [wbdfun, kpfun] = GUI_FreeRotate(ax, enable, coordsys, box)
%
%   Input:
%    ax:       Column vector of handles to the axes on which this function
%              will take effect. These axes should all be on the same
%              figure. Optional, defaults to gca.
%    enable:   Logical indicating whether rotation should be turned on or
%              off. Optional, defaults to false.
%    coordsys: String indicating to which axis the rotation should be
%              constrained. If set to 'none', rotation will not be
%              constrained. If set to 'x', 'y' or 'z', the corresponding
%              axis will be kept vertical. Optional, defaults to 'none'.
%    box:      Logical indicating whether or not plot box rotation should
%              be used instead of continuous rotation. Optional, defaults
%              to false.
%
%   Output:
%    wbdfun: Cell array containing a handle to a locally defined callback
%            function and a set of arguments. This variable can be used as
%            a figure's WindowButtonDownFcn property.
%    kpfun:  Cell array containing a handle to a locally defined callback
%            function and a set of arguments. This variable can be used as
%            a figure's KeyPressFcn property.
%
%   Effect: This function will enable the user to interactively rotate a
%   3-D plot. It also possible to constrain the rotation so the x-, y- or
%   z-axis remains vertical. During constrained rotation, the user may
%   orbit the camera with the left mouse button, zoom with the right mouse
%   button and pan with the middle mouse button. The Page Up and Page Down
%   keys can also be used to zoom in and out. In unconstrained rotation,
%   the right button is used to roll the camera instead of zooming. The
%   camera properties of the axes referenced by ax will be linked together.
%   In order to keep more control over the initialisation and cleanup
%   process, it is possible to set the outputs of GUI_FreeRotate(ax, false)
%   as a figure's WindowButtonDownFcn and KeyPressFcn. This may be useful
%   when writing a toolbar.
%
%   Dependencies: GUI_FreeRotate.m (recursive)
%                 GUI_RunCallback.m
%                 GUI_CameraCenterPixel.m
%
%   Known parents: GUI_FreeRotate.m (recursive)
%                  GUI_VisualisationUI.m

%Created on 17/05/2006 by Ward Bartels.
%WB, 18/05/2006: Added constrained rotation.
%WB, 19/05/2006: Camera can now be lined up with axis.
%WB, 23/05/2006: Added camera rolling.
%WB, 24/05/2006: Added camera panning.
%WB, 26/05/2006: Added keyboard camera zooming.
%WB, 02/06/2006: Added plot box rotation.
%WB, 09/06/2006: Added crosshair to panning.
%WB, 27/06/2006: Added rmb camera zooming.
%WB, 05/07/2006: Added axes linking.
%WB, 25/09/2006: Split off GUI_RunCallback.m.
%WB, 17/07/2007: Improved rotation motion scaling.
%Stabile, fully functional.



%---------------%
% Main function %
%---------------%

%Get handle of the figure to which the axes belongs
if nargin<1
    ax = gca;
end
fig = ancestor(ax(1), 'figure');

%Convert coordsys into argument for camorbit; line up cameras
%<LineUpCamera> <LinkAxes>
if nargin>=3 && ismember(coordsys, {'x' 'y' 'z'})
    LineUpCamera(ax, coordsys);
    coordsys = {'data' coordsys};
else
    LinkAxes(ax);
    coordsys = {'camera'};
end

%Assemble WindowButtonDownFcn and KeyPressFcn properties
%<PlotBoxButtonDown> <ButtonDown> <KeyPress>
if nargin>=4 && box
    wbdfun = {@PlotBoxButtonDown ax coordsys};
else
    wbdfun = {@ButtonDown ax coordsys};
end
kpfun = get(fig, 'KeyPressFcn');
if ~iscell(kpfun)
    kpfun = {kpfun};
end
kpfun = [{@KeyPress ax} kpfun(:).'];

%Determine if rotation should be enabled or disabled
if nargin>=2 && enable
    
    %If free rotation is already active, turn it off <<GUI_FreeRotate.m>>
    if isappdata(fig, mfilename)
        GUI_FreeRotate(ax, false, 'none');
    end
    
    %Turn off other interactive modes and save state; keep uicontrols
    ad.state = uiclearmode(fig, 'docontext', mfilename, ax, false);
    uirestore(ad.state, 'uicontrols');
    setappdata(fig, mfilename, ad);
    
    %Set pointer
    setptr(fig, 'rotate');
    
    %Set window button press and key press callback
    set(fig, 'WindowButtonDownFcn', wbdfun);
    set(fig, 'KeyPressFcn', kpfun);
    
else
    
    %Get application data
    ad = getappdata(fig, mfilename);
    
    %Reset interactive properties and remove application data <ButtonUp>
    if ~isempty(ad)
        if isfield(ad, 'box_ax')
            delete(ad.box_ax(ishandle(ad.box_ax)));
        end
        if isfield(ad, 'cross')
            delete(ad.cross(ishandle(ad.cross)));
        end
        if isfield(ad, 'state')
            uirestore(ad.state, 'nouicontrols');
        else
            ButtonUp(fig, [], ax);
        end
        rmappdata(fig, mfilename);
    end
    
    %Remove button press and key press callbacks if they were set
    wbdget = get(fig, 'WindowButtonDownFcn');
    kpget = get(fig, 'KeyPressFcn');
    if ~isempty(wbdget) && iscell(wbdget) && ...
       (isequal(wbdget{1}, @ButtonDown) || isequal(wbdget{1}, @PlotBoxButtonDown))
        set(fig, 'WindowButtonDownFcn', '');
    end
    if ~isempty(kpget) && iscell(kpget) && isequal(kpget{1}, @KeyPress)
        if numel(kpget)>3
            kpset = kpget(3:end);
        elseif numel(kpget)==3
            kpset = kpget{3};
        else
            kpset = '';
        end
        set(fig, 'KeyPressFcn', kpset);
    end
end



%-----------%
% Callbacks %
%-----------%


%Callback, executed when a key is pressed in the figure
function KeyPress(src, evt, ax, varargin)

%Validate ax <FindAxes>
ax = FindAxes(src, ax);
if isempty(ax), return; end

%Catch (shift-)Page Up and (shift-)Page Down (the shift modifier is used in
%conjunction with the AutoHotkey program, which maps mouse wheel rotation
%to keys)
if isempty(evt.Modifier) || isequal(evt.Modifier, {'shift'})
    if strcmp(evt.Key, 'pageup')       %Page Up - zoom in
        arrayfun(@(x) camzoom(x, 1.2), ax);
        return
    elseif strcmp(evt.Key, 'pagedown') %Page Down - zoom out
        arrayfun(@(x) camzoom(x, 0.8), ax);
        return
    end
end

%Pass on the callback <<GUI_RunCallback.m>>
GUI_RunCallback(src, evt, varargin{:});



%Callback, executed when a mouse button is pressed down in plot box mode
function PlotBoxButtonDown(fig, evt, ax, coordsys)

%Get application data and validate ax <FindAxes>
ax = FindAxes(fig, ax);
if isempty(ax), return; end
ad = getappdata(fig, mfilename);

%Call ButtonUp if another button is still pressed <ButtonUp>
if isfield(ad, 'mode') && ~strcmpi(ad.mode, 'none')
    ButtonUp(fig, evt, ax);
    ad = getappdata(fig, mfilename);
end

%Loop over all axes
for ind = 1:numel(ax)
    
    %Create new axes for outline box
    ad.box_ax(ind) = axes('Parent', get(ax(ind), 'Parent'), 'Visible', 'off', ...
                          'HandleVisibility', 'off', 'Drawmode', 'fast');
    
    %Assemble line coordinates
    if numel(coordsys)<2
        lcoords = [0 0 1 1 0 0 0 1 1 0 NaN 0 0 NaN 1 1 NaN 1 1; ...
                   0 0 0 0 0 1 1 1 1 1 NaN 1 0 NaN 1 0 NaN 1 0; ...
                   1 0 0 1 1 1 0 0 1 1 NaN 0 0 NaN 1 1 NaN 0 0];
    else
        lcoords = [0 0 1 1 0 0 0 1 1 0 0 1 1 0 NaN 1 1; ...
                   0 0 0 0 0 1 1 0 1 0 1 1 1 1 NaN 1 0; ...
                   1 0 0 1 1 1 0 0 0 0 0 0 1 1 NaN 1 1];
        jnd = coordsys{2}-'x'+1; %1 for x, 2 for y and 3 for z
        lcoords([jnd 3],:) = lcoords([3 jnd],:);
    end
    
    %Scale line coordinates so outline box coincides with plot box
    lim(1,:) = get(ax(ind),'XLim');
    lim(2,:) = get(ax(ind),'YLim');
    lim(3,:) = get(ax(ind),'ZLim');
    lcoords = lcoords.*(diff(lim, 1, 2)*ones(1, size(lcoords, 2)))+lim(:,1)*ones(1, size(lcoords, 2));
    
    %Create line object
    ad.box_line(ind) = line(lcoords(1,:), lcoords(2,:), lcoords(3,:), ...
                            'Parent', ad.box_ax(ind), 'Erasemode', 'xor', ...
                            'Visible', 'off', 'HandleVisibility', 'off', 'Clipping', 'off');
end

%Assume camera to be above horizon (will be changed in *ButtonMotion)
ad.box_inverted = false;

%Copy view to ad.box_ax <FollowView>
FollowView(ax, ad.box_ax);
drawnow expose;

%Store application data
ad.orig_ax = ax;
setappdata(fig, mfilename, ad);

%Call ButtonDown <ButtonDown>
ButtonDown(fig, evt, ad.box_ax, coordsys);

%Set WindowButtonUpFcn <ButtonUp>
set(fig, 'WindowButtonUpFcn', {@ButtonUp, ax});



%Callback, executed when a mouse button is pressed down
function ButtonDown(fig, evt, ax, coordsys)

%Set figure units to pixels and get current point
set(fig, 'Units', 'pixels');
pt = get(fig, 'CurrentPoint');

%Get application data and validate ax <FindAxes>
ax = FindAxes(fig, ax);
if isempty(ax), return; end
ad = getappdata(fig, mfilename);

%Call ButtonUp if another button is still pressed <ButtonUp>
if isfield(ad, 'mode') && ~strcmpi(ad.mode, 'none')
    ButtonUp(fig, evt, ax);
    ad = getappdata(fig, mfilename);
end

%Differentiate between selection types
switch get(fig, 'SelectionType')
    case 'normal' %Left button - rotate <RotateButtonMotion>
        ad.mode = 'rotate';
        prop = {@RotateButtonMotion ax coordsys};
    case 'alt'    %Right button - zoom or roll camera
        if numel(coordsys)>=2 %Constrained rotation - zoom <ZoomButtonMotion>
            ad.mode = 'zoom';
            prop = {@ZoomButtonMotion ax coordsys};
        else                  %Unconstrained rotation - roll <RollButtonMotion> <ClickedAxes>
            ad.mode = 'roll'; 
            prop = {@RollButtonMotion ax ClickedAxes(fig, ax, pt) coordsys};
        end
    case 'extend' %Middle button (or left & right) - pan camera <PanButtonMotion>
        ad.mode = 'pan';
        for ind = 1:numel(ax)
            ct = camtarget(ax(ind));
            ad.cross(ind) = line([get(ax(ind), 'XLim') NaN ct(1) ct(1) NaN ct(1) ct(1)], ...
                                 [ct(2) ct(2) NaN get(ax(ind), 'YLim') NaN ct(2) ct(2)], ...
                                 [ct(3) ct(3) NaN ct(3) ct(3) NaN get(ax(ind), 'ZLim')], ...
                                 'Parent', ax(ind), 'Color', [0.2 0 0.5]); %3D cross indicating camera target
        end
        set(ax, 'CameraViewAngleMode', 'manual');
        prop = {@PanButtonMotion ax coordsys};
        setptr(fig, 'closedhand');
    case 'open'   %Double click - reset pan and zoom <LinkAxes>
        set(ax(1), 'CameraTargetMode', 'auto');
        set(ax(1), 'CameraPositionMode', 'auto');
        set(ax(1), 'CameraViewAngleMode', 'auto');
        LinkAxes(ax);
        return
    otherwise     %Other button - do nothing (fail-safe)
        return
end

%Add current point to application data
ad.pt = pt;
setappdata(fig, mfilename, ad);

%Add callbacks to figure <ButtonUp>
set(fig, 'WindowButtonMotionFcn', prop);
set(fig, 'WindowButtonUpFcn', {@ButtonUp ax});



%Callback, executed when the mouse cursor is moved within the figure with
%the left button down
function RotateButtonMotion(fig, evt, ax, coordsys)

%Get correct axes handle and motion in pixels <InitButtonMotion>
[ax, startpt, endpt] = InitButtonMotion(fig, ax, coordsys);
if isempty(ax), return; end

%Determine how much the cursor has moved
xy = startpt-endpt;

%Line up camera and invert first element of xy if necessary <LineUpCamera>
if length(coordsys)>=2 && LineUpCamera(ax, coordsys{2})
    xy(1) = -xy(1);
end

%Get axes size in pixels
pos = hgconvertunits(fig, get(ax(1), 'Position'), get(ax(1), 'Units'), 'pixels', get(ax(1), 'Parent'));
% units = get(ax(1), 'Units');
% set(ax(1), 'Units', 'pixels');
% pos = get(ax(1), 'Position');
% set(ax(1), 'Units', units);
pix = min(pos(3), pos(4));

%Scale movement so cursor "holds on" to outer boundary of plot box
siz = min(abs(diff([xlim(ax(1)); ylim(ax(1)); zlim(ax(1))], [], 2)))/2; %plot box minimum extents
dis = norm((camtarget(ax(1))-campos(ax(1)))./daspect(ax(1)));           %distance from camera target to camera
fov = 2*dis*tand(camva(ax(1))/2);                                       %field of view in axes coords
% xy = atand(xy*fov/pix/siz); %exact "rolling" motion
xy = xy*fov/pix/siz*180/pi;   %slightly accelerated motion

%Orbit the camera and redraw
arrayfun(@(x) camorbit(x, xy(1), xy(2), coordsys{:}), ax);
%drawnow expose; %This is used in cameratoolbar.m and rotate3d.m to fix a bug



%Callback, executed when the mouse cursor is moved within the figure with
%the right button down and constrained rotation active
function ZoomButtonMotion(fig, evt, ax, coordsys)

%Get correct axes handle and motion in pixels <InitButtonMotion>
[ax, startpt, endpt, ad] = InitButtonMotion(fig, ax, coordsys);
if isempty(ax), return; end

%Get size of figure window
figsize = get(fig, 'Position');
figsize = figsize(3:4);

%Determine how much the cursor has moved
xy = endpt-startpt;

%Normalize xy to figsize and cap at -1 and 1
xy = xy./figsize;
xy = max(min(xy, 1), -1);

%Check if y-motion is largest
if abs(xy(1))<abs(xy(2))
    
    %Set the pointer to zoom in or out
    if xy(2)>0
        setptr(fig, 'glassplus');
    else
        setptr(fig, 'glassminus');
    end
    
    %Zoom the camera
    arrayfun(@(x) camzoom(x, 100^xy(2)), ax);
end



%Callback, executed when the mouse cursor is moved within the figure with
%the right button down and free rotation active
function RollButtonMotion(fig, evt, ax, ax_click, coordsys)

%Get correct axes handle and motion in pixels <InitButtonMotion>
[ax, startpt, endpt, ad] = InitButtonMotion(fig, ax, coordsys);
if isempty(ax), return; end

%Get the camera target in pixel coordinates <<GUI_CameraCenterPixel.m>>
center = GUI_CameraCenterPixel(ax_click);

%Calculate rotation in degrees
v1 = startpt-center;
v2 = endpt-center;
theta = real(acos(v1*v2.'/norm(v1, 2)/norm(v2, 2))*180/pi);

%Check if rotation is clockwise
if v2(2)*v1(1)<v2(1)*v1(2)
    theta = -theta;
end

%Roll the camera
arrayfun(@(x) camroll(x, theta), ax);
%drawnow expose; %This is used in cameratoolbar.m and rotate3d.m to fix a bug



%Callback, executed when the mouse cursor is moved within the figure with
%the right button down and constrained rotation active
function PanButtonMotion(fig, evt, ax, coordsys)

%Get correct axes handle and motion in pixels <InitButtonMotion>
[ax, startpt, endpt, ad] = InitButtonMotion(fig, ax, coordsys);
if isempty(ax), return; end

%Dolly the camera
arrayfun(@(x) camdolly(x, startpt(1)-endpt(1), startpt(2)-endpt(2), 0, 'movetarget', 'pixels'), ax);

%Update cross position
if isfield(ad, 'cross') && all(ishandle(ad.cross))
    for ind = 1:numel(ax)
        ct = camtarget(ax(ind));
        set(ad.cross(ind), ...
            'XData', [get(ax(ind), 'XLim') NaN ct(1) ct(1) NaN ct(1) ct(1)], ...
            'YData', [ct(2) ct(2) NaN get(ax(ind), 'YLim') NaN ct(2) ct(2)], ...
            'ZData', [ct(3) ct(3) NaN ct(3) ct(3) NaN get(ax(ind), 'ZLim')]);
    end
end
%drawnow expose; %This is used in cameratoolbar.m and rotate3d.m to fix a bug



%Callback, executed when the mouse button is released
function ButtonUp(fig, evt, ax)

%Remove callbacks (up and motion)
set(fig, 'WindowButtonMotionFcn', '');
set(fig, 'WindowButtonUpFcn', '');

%Get application data and validate ax <FindAxes>
ad = getappdata(fig, mfilename);
ax = FindAxes(fig, ax);

%Fail-safe
if isempty(ax)
    return
end

%If mode is 'pan', reset pointer
if any(strcmp(ad.mode, {'pan' 'zoom'}))
    if isfield(ad, 'cross') 
        delete(ad.cross(ishandle(ad.cross)));
        ad = rmfield(ad, 'cross');
    end
    setptr(fig, 'rotate');
end

%Remove outline box if it is present; copy view to ax <FollowView>
if isfield(ad, 'box_ax')
    if all(ishandle(ad.box_ax))
        FollowView(ad.box_ax, ax);
        delete(ad.box_ax);
    end
    ad = rmfield(ad, {'orig_ax' 'box_ax' 'box_line' 'box_inverted'});
end

%Reset mode
ad.mode = 'none';
setappdata(fig, mfilename, ad);



%-------------------%
% Utility functions %
%-------------------%


%Utility function, initialises a ButtonMotion function
function [ax, startpt, endpt, ad] = InitButtonMotion(fig, ax, coordsys)

%If ax no longer exists, find a new axes <FindAxes>
ax = FindAxes(fig, ax);

%Get application data
ad = getappdata(fig, mfilename);

%Fail-safe (should make caller function return)
if isempty(ad) || isempty(ax)
    startpt = [0 0];
    endpt = [0 0];
    ax = [];
    return
end

%Check if outline box is present
if isfield(ad, 'box_line') && all(ishandle(ad.box_line))
    
    %Make outline box visible
    set(ad.box_line, 'Visible', 'on');
    
    %Check if rotation is constrained
    if numel(coordsys)>=2
        
        %Get axis index
        ind = coordsys{2}-'x'+1; %1 for x, 2 for y and 3 for z
        
        %Check if camera is below horizon
        camvect = campos(ax(1))-camtarget(ax(1));
        cam_inverted = camvect(ind)<0;
        
        %Invert box if necessary
        if xor(cam_inverted, ad.box_inverted)
            prop = [upper(coordsys{2}) 'Data'];
            jnd = [1 4:6 13 14 16 17; 2 3 7:12];
            for knd = 1:numel(ad.box_line)
                coord = get(ad.box_line(knd), prop);
                coord([jnd(1,:) jnd(2,:)]) = coord([jnd(2,:) jnd(1,:)]);
                set(ad.box_line(knd), prop, coord);
            end
            ad.box_inverted = ~ad.box_inverted;
        end
    end
end

%Determine start and end point of motion, update current point
set(fig, 'Units', 'pixels');
endpt = get(fig, 'CurrentPoint');
startpt = ad.pt;
ad.pt = endpt;
setappdata(fig, mfilename, ad);



%Utility function, finds a new axes if ax no longer exists
function ax = FindAxes(fig, ax)

%Eliminate invalid axes
ax(~ishandle(ax)) = [];

%If no valid axes exists, find a new axes
if isempty(ax)
    ax = findobj(fig, 'Type', 'axes');
end



%Utility function, finds out which axes was clicked
function ax = ClickedAxes(fig, ax, startpt)

%Get axes positions in pixel coordinates
units = get(ax, {'Units'});
set(ax, 'Units', 'pixels');
pos = get(ax, {'Position'});
set(ax, {'Units'}, units);
pos = vertcat(pos{:});

%Check which axes contain the clicked point
ind = all([startpt(1)>=pos(:,1) startpt(2)>=pos(:,2) ...
           startpt(1)<=pos(:,1)+pos(:,3) startpt(2)<=pos(:,2)+pos(:,4)], 2);
if sum(ind)==0
    ind = true(size(ind));
end
ax = ax(ind);

%If multiple axes remain, check which center is closest to clicked point
if numel(ax)>1
    pos = pos(ind,:);
    center = pos(:,1:2)+pos(:,3:4)/2;
    dist = (center(:,1)-startpt(1)).^2+center(:,2)-startpt(2).^2;
    [ignoble, ind] = min(dist);
    ax = ax(ind);
end



%Utility function, lines up the camera with the selected axis
function inverted = LineUpCamera(ax, coordsys)

%Create axis vector corresponding to coordsys
ind = coordsys-'x'+1; %1 for x, 2 for y and 3 for z
axisvect = [0 0 0];
axisvect(ind) = 1;

%If view is upside-down, change axisvect direction
upvect = camup(ax(1));
inverted = upvect(ind)<0;
if inverted
    axisvect = -axisvect;
end

%Set the camera up vector if axisvect is not parallel with view direction
if any(cross(axisvect, campos(ax(1))-camtarget(ax(1))))
    camup(ax(1), axisvect);
end

%Copy view from ax(1) to all other axes <LinkAxes>
LinkAxes(ax);



%Utility function, copies camera properties from first axes to others
function LinkAxes(ax)

%Copy view from ax(1) to all other axes
properties = {'CameraPosition' 'CameraTarget' 'CameraTargetMode' ...
              'CameraViewAngle' 'CameraViewAngleMode' 'CameraUpVector' ...
              'CameraUpVectorMode'};
set(ax(2:end), properties, get(ax(1), properties));



%Utility function, copies view properties from one axes to another
function FollowView(ax_src, ax_dst)

%List of relevant properties
properties = {'XLim' 'YLim' 'ZLim' 'DataAspectRatio' 'DataAspectRatioMode' ...
              'PlotBoxAspectRatio' 'PlotBoxAspectRatioMode' ...
              'CameraPosition' 'CameraTarget' 'CameraTargetMode' ...
              'CameraViewAngle' 'CameraViewAngleMode' 'CameraUpVector' ...
              'CameraUpVectorMode' 'Projection' 'Units' 'Position'};

%Copy properties
set(ax_dst, properties, get(ax_src, properties));