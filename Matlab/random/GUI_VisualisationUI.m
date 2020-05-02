function [tb, fun] = GUI_VisualisationUI(fig, removeui, ax, toolbar, keys)
%GUI_VisualisationUI  Add user interface functions to visualisations.
%wb20060705
%
%   Syntax:
%    [tb, fun] = GUI_VisualisationUI(fig, removeui, ax, toolbar, keys)
%
%   Input:
%    fig:      Figure handle.
%    removeui: Logical indicating whether or not the original toolbar and
%              menu bar should be removed.
%    ax:       Column vector of handles to the axes on which this function
%              will take effect. If ax is empty, the axes under fig will be
%              used. Optional, only required when either toolbar or keys is
%              set to true.
%    toolbar:  Logical indicating whether or not the visualisation toolbar
%              should be added to the figure. Optional, defaults to false.
%    keys:     Logical indicating whether or not shortcut keys should be
%              activated. Optional, defaults to false.
%
%   Output:
%    tb:  Column vector containing handles to the toolbar and its children.
%         The first element is a handle to the toolbar, the next elements
%         are handles to the toolbar's children. If toolbar is set to
%         false, tb will be an empty matrix.
%    fun: Cell array containing a handle to GUI_VisualisationKeyPress.m and
%         a set of arguments. This variable can be used as a figure's
%         KeyPressFcn property. If keys is set to false, fun will be an
%         empty cell array.
%
%   Effect: This function can remove the original toolbar and menu bar from
%   a figure, add a standard visualisation toolbar and activate shortcut
%   keys, depending on the value of removeui, toolbar and keys. If a
%   toolbar with tag 'VisToolbar' already exists in fig, the buttons will
%   be added to this toolbar; otherwise, a new toolbar will be created with
%   this tag. The toolbar allows the user to orbit the camera in two ways
%   (continuous orbit and plot box orbit), constrain orbiting to keep the
%   Z-axis vertical, line up the lights with the camera, toggle perspective
%   and toggle visibility of the plot box. The shortcut keys are the ones
%   set in GUI_VisualisationKeyPress.m.
%
%   Dependencies: GUI_LoadPNGButtons.m
%                 GUI_VisualisationKeyPress.m
%                 GUI_FreeRotate.m
%                 Orbit.png
%                 View.png
%                 Free rotate.png
%                 Light.png
%                 Perspective.png
%                 Box.png
%
%   Known parents: GUI_SimulateSurgery.m

%Created on 24/04/2006 by Ward Bartels.
%WB, 27/04/2006: Added cursor to free rotation.
%WB, 28/04/2006: Improved interaction with other interactive modes.
%WB, 03/05/2006: Deleting the axes will no longer cause errors; axes handle
%                can now be determined automatically.
%WB, 05/05/2006: Other interactive modes will now be turned off first;
%                existing toolbar is now reused thanks to Tag properties.
%WB, 02/06/2006: GUI_FreeRotate.m is now used instead of cameratoolbar and
%                rotate3d.
%WB, 05/07/2006: Added axes linking.
%Stabile, fully functional.



%---------------%
% Main function %
%---------------%


%Turn off other interactive modes
scribeclearmode(fig);

%Remove original UI elements if requested
if removeui
    set(fig, 'ToolBar', 'none', 'MenuBar', 'none');
end

%Add shortcut keys if requested <<GUI_VisualisationKeyPress.m>>
if nargin>=5 && keys
    fun = {@GUI_VisualisationKeyPress ax};
    set(fig, 'KeyPressFcn', fun);
else
    fun = {};
end

%Add toolbar if requested
if nargin>=4 && toolbar
    
    %Load toolbar button icons <<GUI_LoadPNGButtons.m>>
    [imorb, imview, imzcon, imbox, imucon, imlight, improj] = ...
    GUI_LoadPNGButtons([], 'Orbit.png', 'View.png', 'Z constrained.png', 'Box.png', ...
                           'Unconstrained.png', 'Light.png', 'Perspective.png');
    
	%Store tags
    tags = {'OrbTool'; 'ViewTool'; 'ConTool'; 'LightTool'; 'ProjTool'; 'BoxTool'};
                           
    %Create toolbar, overwrite existing toolbar with tag VisToolbar
    tb = findobj(fig, 'Tag', 'VisToolbar');
    tb(2:end) = [];
    if isempty(tb)
        tb(1) = uitoolbar('Parent', fig, 'Tag', 'VisToolbar');
    else
        tbchild = get(tb, 'Children');
        delete(tbchild(ismember(get(tbchild, 'Tag'), tags)));
    end
    
    %Create buttons <OrbitOn> <OrbitViewOff> <ViewOn> <FRotOn> <FRotOff>
    %               <LightTool> <ProjOn> <ProjOff> <BoxOn> <BoxOff>
    tb(2,1) = uitoggletool('Parent', tb(1), 'CData', imorb, 'Tag', tags{1}, 'TooltipString', 'Continuous orbit', ...
                           'OnCallback', {@OrbitViewOn, ax, fig}, 'OffCallback', {@OrbitViewOff, ax, fig});
	tb(3,1) = uitoggletool('Parent', tb(1), 'CData', imview, 'Tag', tags{2}, 'TooltipString', 'Plot box orbit', ...
                           'OnCallback', {@OrbitViewOn, ax, fig}, 'OffCallback', {@OrbitViewOff, ax, fig});
	tb(4,1) = uipushtool('Parent', tb(1), 'CData', imzcon, 'Tag', tags{3}, 'TooltipString', 'Change rotation constraint', ...
                         'ClickedCallback', @ConTool, 'UserData', {imucon 'none'; imzcon 'z'});
	tb(5,1) = uipushtool('Parent', tb(1), 'CData', imlight, 'Tag', tags{4}, 'TooltipString', 'Align lights with camera', ...
                         'ClickedCallback', {@LightTool, ax, fig}, 'Separator', 'on');
	tb(6,1) = uitoggletool('Parent', tb(1), 'CData', improj, 'Tag', tags{5}, 'TooltipString', 'Toggle perspective', ...
                           'OnCallback', {@ProjOn, ax, fig}, 'OffCallback', {@ProjOff, ax, fig});
	tb(7,1) = uitoggletool('Parent', tb(1), 'CData', imbox, 'Tag', tags{6}, 'TooltipString', 'Toggle plot box', ...
                           'OnCallback', {@BoxOn, ax, fig}, 'OffCallback', {@BoxOff, ax, fig});
	
	%Make sure toggle buttons are turned off on deletion
    set(tb(2:3), {'DeleteFcn'}, get(tb(2:3), {'OffCallback'}));
    
    %Set state of last two toggle tools <FindAxes>
    ax = FindAxes(fig, ax);
    if strcmp(camproj(ax(1)), 'perspective')
        set(tb(6), 'State', 'on');
    end
    if strcmp(get(ax(1), 'Visible'), 'on')
        set(tb(7), 'State', 'on');
    end
    
    %Mark OrbTool as last used and set rotation constraint to Z-axis
    setappdata(tb(1), 'LastViewButton', tb(2));
    setappdata(tb(1), 'RotationConstraint', 'z');
    
else
    
    %Pass empty matrix to output
    tb = [];
end



%-----------%
% Callbacks %
%-----------%


%Callback, executed when the Orbit or View tool is activated
function OrbitViewOn(src, evt, ax, fig)

%Validate axes <FindAxes>
ax = FindAxes(fig, ax);
if isempty(ax), return; end

%Check which tool was activated
box = strcmpi(get(src, 'Tag'), 'ViewTool');

%Check rotation constraint setting
coordsys = getappdata(get(src, 'Parent'), 'RotationConstraint');

%Turn off other toolbar buttons
scribeclearmode(fig);

%Store callbacks
state = uiclearmode(fig, 'docontext', @set, src, 'State', 'off');
uirestore(state, 'uicontrols');
setappdata(fig, 'GUI_VisualisationUI', state);

%Enable free rotation <<GUI_FreeRotate.m>>
[wbdfun, kpfun] = GUI_FreeRotate(ax, false, coordsys, box);
set(fig, 'WindowButtonDownFcn', wbdfun);
set(fig, 'KeyPressFcn', kpfun);

%Set pointer
setptr(fig, 'rotate');

%Mark button as last used
setappdata(get(src, 'Parent'), 'LastViewButton', src);



%Callback, executed when the Orbit or View tool is deactivated
function OrbitViewOff(src, evt, ax, fig)

%Validate axes <FindAxes>
ax = FindAxes(fig, ax);
if isempty(ax), return; end

%Disable free rotation <<GUI_FreeRotate.m>>
GUI_FreeRotate(ax, false, 'none', false);

%Restore interactive modes
if isappdata(fig, 'GUI_VisualisationUI')
    uirestore(getappdata(fig, 'GUI_VisualisationUI'), 'nouicontrols');
    rmappdata(fig, 'GUI_VisualisationUI')
end

%Remove ScribeClearModeCallback if it exists
if isappdata(fig, 'ScribeClearModeCallback')
    rmappdata(fig, 'ScribeClearModeCallback')
end



%Callback, executed when the user clicks on the Constraint toolbar button
function ConTool(src, evt)

%Get UserData (images and corresponding constraint string)
ud = get(src, 'UserData');

%Check which row of ud is currently selected; increment
ind = strmatch(getappdata(get(src, 'Parent'), 'RotationConstraint'), ud(:,2));
if isempty(ind)
    ind = 1;
end
ind = mod(ind, size(ud, 1))+1;

%Set application data and CData accordingly
setappdata(get(src, 'Parent'), 'RotationConstraint', ud{ind,2});
set(src, 'CData', ud{ind,1});



%Callback, executed when the user clicks on the Light toolbar button
function LightTool(src, evt, ax, fig)

%Validate axes <FindAxes>
ax = FindAxes(fig, ax);
if isempty(ax), return; end

%Loop over all axes
for ind = 1:numel(ax)
    
    %Get handles to light objects
    li = get(ax(ind), 'Children');
    li = li(strcmp(get(li, 'Type'), 'light'));
    
    %Move lights
    if length(li)>=1
        camlight(li(1), 'headlight', 'infinite');
    end
    if length(li)>=2
        camlight(li(2), 0, 180, 'infinite');
    end
end



%Callback, executed when the Perspective tool is activated
function ProjOn(src, evt, ax, fig)

%Validate axes <FindAxes>
ax = FindAxes(fig, ax);
if isempty(ax), return; end

%Set camera projection to perspective
set(ax, 'Projection', 'perspective');



%Callback, executed when the Perspective tool is deactivated
function ProjOff(src, evt, ax, fig)

%Validate axes <FindAxes>
ax = FindAxes(fig, ax);
if isempty(ax), return; end

%Set camera projection to orthographic
set(ax, 'Projection', 'orthographic');



%Callback, executed when the Plot box tool is activated
function BoxOn(src, evt, ax, fig)

%Validate axes <FindAxes>
ax = FindAxes(fig, ax);
if isempty(ax), return; end

%Set ax visible
set(ax, 'Visible', 'on');



%Callback, executed when the Plot box tool is deactivated
function BoxOff(src, evt, ax, fig)

%Validate axes <FindAxes>
ax = FindAxes(fig, ax);
if isempty(ax), return; end

%Set ax invisible
set(ax, 'Visible', 'off');



%-------------------%
% Utility functions %
%-------------------%


%Utility function, finds the first axes under a figure
function ax = FindAxes(fig, ax)

%Fail-safe, in case toolbar is moved to another figure and fig is deleted
if ~ishandle(fig)
    ax = [];
    return
end

%Eliminate invalid axes
ax(~ishandle(ax)) = [];

%If no valid axes exists, find a new axes
if isempty(ax)
    ax = findobj(fig, 'Type', 'axes');
end