function varargout = GUI_MeshExplorer(varargin)
%GUI_MeshExplorer  Display a mesh interactively.
%wb20060315
%
%   Syntax:
%    fig = GUI_MeshExplorer(ax)
%
%   Input:
%    ax: Handle to the axes containing the objects to be displayed.
%
%   Output:
%    fig: Handle to the figure window created by GUI_MeshExplorer.
%
%   Effect: This function will copy the objects contained within ax to a
%   new window. The user can interactively manipulate the camera
%   orientation, distance and view angle using buttons displayed in the
%   window. Perspective can be switched on and off. It is also possible to
%   switch to a crossed-eye or a parallel-eye stereoscopic view. When
%   steroscopy is active, the intraocular distance coefficient (IDC) can be
%   changed. The distance between the two cameras is calculated as D/IDC,
%   D being the distance from one of the cameras to its target. The
%   following keys can be used:
%      - Up arrow:    Tilt mesh upward
%      - Down arrow:  Tilt mesh downward
%      - Left arrow:  Rotate mesh towards the left
%      - Right arrow: Rotate mesh towards the right
%      - 0:           Zoom out (increase view angle)
%      - 1:           Zoom in (decrease view angle)
%      - 4:           Walk out (increase camera distance)
%      - 7:           Walk in (decrease camera distance)
%      - Enter:       Reset camera distance and view angle
%      - Space:       Show or hide the control panel
%   As soon as the window is created, GUI_MeshExplorer exits and returns a
%   handle to the window. Use uiwait(handle) to suspend execution until the
%   window is closed.
%
%   Dependencies: GUI_MeshExplorer.fig
%
%   Known parents: GUI_VisualisationKeyPress.m

%Created on 10/03/2006 by Ward Bartels with GUIDE v2.5.
%WB, 15/03/2006: Added fail-safes to prevent problems when non-numerical
%                text is entered into IDC or when the window is resized to
%                a height smaller than the control panel.
%Stabile, fully functional.



%----------------%
% Initialisation %
%----------------%


%Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_MeshExplorer_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_MeshExplorer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
%End initialization code - DO NOT EDIT



%Executed just before GUI_MeshExplorer is made visible
function GUI_MeshExplorer_OpeningFcn(hObject, eventdata, handles, varargin)
% hObject    handle to Window
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_MeshExplorer (see VARARGIN)

%Choose default command line output for GUI_MeshExplorer
handles.output = hObject;

%Get handle of original axes
if nargin>3
    ax_orig = varargin{1};
else
    error('GUI_MeshExplorer: Not enough input arguments.');
end

%Get handles to new axes
ax_new = [handles.Left_Axis; handles.Center_Axis; handles.Right_Axis];

%Copy data from original axes to new axes
obj_orig = get(ax_orig, 'Children');
copyobj(obj_orig, ax_new(1));
copyobj(obj_orig, ax_new(2));
copyobj(obj_orig, ax_new(3));

%Set objects on left and right axes invisible
set(get(ax_new(1), 'Children'), 'Visible', 'off');
set(get(ax_new(3), 'Children'), 'Visible', 'off');

%Set axes properties
axis(ax_new, 'tight');
set(ax_new, 'CameraPositionMode', 'auto', 'CameraTargetMode', 'auto', 'CameraUpVectorMode', 'auto', 'CameraViewAngleMode', 'auto');

%Copy camera orientation from original axes to new axes
[az, el] = view(ax_orig);
view(ax_new(2), [az el]);

%Save camera orientation
set(hObject, 'UserData', struct('az', az, 'el', el));

%Update handles structure
guidata(hObject, handles);



%Outputs from this function are returned to the command line
function varargout = GUI_MeshExplorer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get default command line output from handles structure
varargout{1} = handles.output;



%-----------%
% Callbacks %
%-----------%


%Callback, executed during object creation, after setting all properties
function IDC_CreateFcn(hObject, eventdata, handles)
%hObject    handle to IDC (see GCBO)
%eventdata  reserved - to be defined in a future version of MATLAB
%handles    empty - handles not created until after all CreateFcns called

%Set background to white
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%Callback, executed when the value in IDC is changed
function IDC_Callback(hObject, eventdata, handles)
% hObject    handle to IDC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Use stereo callbacks <Stereo_Crossed_Callback>
%                     <Stereo_Parallel_Callback>
Stereo_Crossed_Callback(handles.Stereo_Crossed, eventdata, handles);
Stereo_Parallel_Callback(handles.Stereo_Parallel, eventdata, handles);



%Callback, executed on button press in Perspective
function Perspective_Callback(hObject, eventdata, handles)
% hObject    handle to Perspective (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Set projection
if get(hObject,'Value')
    for ax = [handles.Left_Axis handles.Center_Axis handles.Right_Axis]
        camproj(ax, 'perspective');
    end
else
    for ax = [handles.Left_Axis handles.Center_Axis handles.Right_Axis]
        camproj(ax, 'orthographic');
    end
end



%Callback, executed on button press in Up
function Up_Callback(hObject, eventdata, handles)
% hObject    handle to Up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Decrease elevation for visible axes <VisibleAxes>
for ax = VisibleAxes(handles).'
    camorbit(ax, 0, -5);
end



%Callback, executed on button press in Down
function Down_Callback(hObject, eventdata, handles)
% hObject    handle to Down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Increase elevation for visible axes <VisibleAxes>
for ax = VisibleAxes(handles).'
    camorbit(ax, 0, 5);
end



%Callback, executed on button press in Center
function Center_Callback(hObject, eventdata, handles)
% hObject    handle to Center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Read saved state
st = get(handles.Window, 'UserData');

%Get handles of visible axes <VisibleAxes>
ax = VisibleAxes(handles);

%Distinguish between stereoscopic and non-stereoscopic
if length(ax)==1 %Non-stereo
    
    %Return camera orientation to saved state
    view(ax, [st.az st.el]);
    
else             %Stereo
    
    %Get azimuth and elevation of both axes (azimuth difference is needed)
    [az1, el1] = view(ax(1));
    [az2, el2] = view(ax(2));
    
    %Return camera orientation to saved state
    view(ax(1), [st.az st.el]);
    view(ax(2), [st.az+az2-az1 st.el]);
end

%Reset camera up vector
set(ax, 'CameraUpVectorMode', 'auto');



%Callback, executed on button press in Right
function Right_Callback(hObject, eventdata, handles)
% hObject    handle to Right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Decrease azimuth for visible axes <VisibleAxes>
for ax = VisibleAxes(handles).'
    camorbit(ax, -5, 0);
end



%Callback, executed on button press in Left
function Left_Callback(hObject, eventdata, handles)
% hObject    handle to Left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Increase azimuth for visible axes <VisibleAxes>
for ax = VisibleAxes(handles).'
    camorbit(ax, 5, 0);
end



%Callback, executed on button press in Zoom_In
function Zoom_In_Callback(hObject, eventdata, handles)
% hObject    handle to Zoom_In (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Zoom in visible axes <VisibleAxes>
for ax = VisibleAxes(handles).'
    camzoom(ax, 1.1);
end



%Callback, executed on button press in Zoom_Reset
function Zoom_Reset_Callback(hObject, eventdata, handles)
% hObject    handle to Zoom_Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Reset zoom for visible axes <VisibleAxes>
set(VisibleAxes(handles), 'CameraViewAngleMode', 'auto');



%Callback, executed on button press in Zoom_Out
function Zoom_Out_Callback(hObject, eventdata, handles)
% hObject    handle to Zoom_Out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Zoom out visible axes <VisibleAxes>
for ax = VisibleAxes(handles).'
    camzoom(ax, 0.9);
end



%Callback, executed on button press in Walk_In
function Walk_In_Callback(hObject, eventdata, handles)
% hObject    handle to Walk_In (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Decrease walk distance for visible axes <VisibleAxes>
for ax = VisibleAxes(handles).'
    
    %Get camera position and target
    pos = campos(ax);
    tar = camtarget(ax);
    
    %Move camera closer
    campos(ax, tar+(pos-tar)*0.9);
end



%Callback, executed on button press in Walk_Reset
function Walk_Reset_Callback(hObject, eventdata, handles)
% hObject    handle to Walk_Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Reset walk distance for visible axes <VisibleAxes>
for ax = VisibleAxes(handles).'
    
    %Save camera orientation
    [az, el] = view(ax);
    
    %Reset camera position and up vector
    set(ax, 'CameraPositionMode', 'auto')
    set(ax, 'CameraUpVectorMode', 'auto');
    
    %Return camera orientation to saved state
    view(ax, [az el]);
end



%Callback, executed on button press in Walk_Out
function Walk_Out_Callback(hObject, eventdata, handles)
% hObject    handle to Walk_Out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Increase walk distance for visible axes <VisibleAxes>
for ax = VisibleAxes(handles).'
    
    %Get camera position and target
    pos = campos(ax);
    tar = camtarget(ax);
    
    %Move camera further
    campos(ax, tar+(pos-tar)*1.1);
end



%Callback, executed on button press in Stereo_Off
function Stereo_Off_Callback(hObject, eventdata, handles)
% hObject    handle to Stereo_Off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Only run the following if Stereo_Off is on
if get(hObject,'Value')
    
    %Get handle of first visible axes <VisibleAxes>
    ax = VisibleAxes(handles);
    ax = ax(1);
    
    %Get handles of left, right and center axes
    ax_left   = handles.Left_Axis;
    ax_center = handles.Center_Axis;
    ax_right  = handles.Right_Axis;
    
    %Copy camera orientation from ax to ax_right
    campos(ax_center, campos(ax));
    camtarget(ax_center, camtarget(ax));
    
    %Copy camera view angle from ax to ax_center
    mode  = camva(ax, 'mode');
    angle = camva(ax);
    camva(ax_center, angle);
    camva(ax_center, mode);
    
    %Set background color to white
    set(handles.Window, 'Color', [224 223 227]/255);
    
    %Set left and right axes invisible, and center axes visible
    set(ax_left, 'Visible', 'off');
    set(get(ax_left, 'Children'), 'Visible', 'off');
    set(ax_right, 'Visible', 'off');
    set(get(ax_right, 'Children'), 'Visible', 'off');
    set(ax_center, 'Visible', 'on');
    set(get(ax_center, 'Children'), 'Visible', 'on');
    
    %Reset camera up vector
    set(ax_center, 'CameraUpVectorMode', 'auto');
end



%Callback, executed on button press in Stereo_Crossed
function Stereo_Crossed_Callback(hObject, eventdata, handles)
% hObject    handle to Stereo_Crossed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get IDC
idc = str2double(get(handles.IDC, 'String'));

%Only run the following if Stereo_Crossed is on and idc is valid
if get(hObject,'Value') && ~isnan(idc)
    
    %Get handle of first visible axes <VisibleAxes>
    ax = VisibleAxes(handles);
    ax = ax(1);
    
    %Get handles of left, right and center axes
    ax_left   = handles.Left_Axis;
    ax_center = handles.Center_Axis;
    ax_right  = handles.Right_Axis;
    
    %Copy camera orientation from ax to ax_right
    campos(ax_right, campos(ax));
    camtarget(ax_right, camtarget(ax));
    
    %Copy camera view angle from ax to ax_left and ax_right
    mode  = camva(ax, 'mode');
    angle = camva(ax);
    camva(ax_left, angle);
    camva(ax_left, mode);
    camva(ax_right, angle);
    camva(ax_right, mode);
    
    %Separate cameras <Stereo>
    Stereo(ax_right, ax_left, idc);
    
    %Set background color to white
    set(handles.Window, 'Color', 'white');
    
    %Set left and right axes visible, and center axes invisible
    set(ax_center, 'Visible', 'off');
    set(get(ax_center, 'Children'), 'Visible', 'off');
    set(ax_left, 'Visible', 'on');
    set(get(ax_left, 'Children'), 'Visible', 'on');
    set(ax_right, 'Visible', 'on');
    set(get(ax_right, 'Children'), 'Visible', 'on');
end



%Callback, executed on button press in Stereo_Parallel
function Stereo_Parallel_Callback(hObject, eventdata, handles)
% hObject    handle to Stereo_Parallel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get IDC
idc = str2double(get(handles.IDC, 'String'));

%Only run the following if Stereo_Parallel is on and idc is valid
if get(hObject,'Value') && ~isnan(idc)
    
    %Get handle of first visible axes <VisibleAxes>
    ax = VisibleAxes(handles);
    ax = ax(1);
    
    %Get handles of left, right and center axes
    ax_left   = handles.Left_Axis;
    ax_center = handles.Center_Axis;
    ax_right  = handles.Right_Axis;
    
    %Copy camera orientation from ax to ax_right
    campos(ax_left, campos(ax));
    camtarget(ax_left, camtarget(ax));
    
    %Copy camera view angle from ax to ax_left and ax_right
    mode  = camva(ax, 'mode');
    angle = camva(ax);
    camva(ax_left, angle);
    camva(ax_left, mode);
    camva(ax_right, angle);
    camva(ax_right, mode);
    
    %Separate cameras <Stereo>
    Stereo(ax_left, ax_right, idc);
    
    %Set background color to white
    set(handles.Window, 'Color', 'white');
    
    %Set left and right axes visible, and center axes invisible
    set(ax_center, 'Visible', 'off');
    set(get(ax_center, 'Children'), 'Visible', 'off');
    set(ax_left, 'Visible', 'on');
    set(get(ax_left, 'Children'), 'Visible', 'on');
    set(ax_right, 'Visible', 'on');
    set(get(ax_right, 'Children'), 'Visible', 'on');
end



%Callback, executed on key press over Window with no controls selected
function Window_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to Window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Determine which key was pressed
switch get(hObject, 'CurrentCharacter')
    case char(30), %Up arrow - decrease elevation <Up_Callback>
        Up_Callback([], [], handles);
    case char(31), %Down arrow - increase elevation <Down_Callback>
        Down_Callback([], [], handles);
    case char(29), %Right arrow - decrease azimuth <Right_Callback>
        Right_Callback([], [], handles);
    case char(28), %Left arrow - increase azimuth <Left_Callback>
        Left_Callback([], [], handles);
    case '1'       %1 - zoom in <Zoom_In_Callback>
        Zoom_In_Callback([], [], handles);
    case '0'       %0 - zoom out <Zoom_Out_Callback>
        Zoom_Out_Callback([], [], handles);
    case char(13)  %Enter - reset zoom & walk <Zoom_Reset_Callback>
                   %                          <Walk_Reset_Callback>
        Zoom_Reset_Callback([], [], handles);
        Walk_Reset_Callback([], [], handles);
    case '7'       %7 - walk in <Walk_In_Callback>
        Walk_In_Callback([], [], handles);
    case '4'%      %4 - walk out <Walk_Out_Callback>
        Walk_Out_Callback([], [], handles);
    case ' '       %Space - show/hide panel <TogglePanel>
        TogglePanel(handles);
end



%Callback, executed on key press over Window with a control selected
function KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to Window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Call Window's KeyPressFcn <Window_KeyPressFcn>
Window_KeyPressFcn(handles.Window, eventdata, handles);



%Callback, executed when Window is resized
function Window_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to Window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get axes handles
ax = [handles.Left_Axis; handles.Center_Axis; handles.Right_Axis];

%Original position for Axis_Left, Axis_Center and Axis_Right (normalized)
positions = [get(ax(1), 'Position'); get(ax(2), 'OuterPosition'); get(ax(3), 'Position')];

%Get window height
pos = get(handles.Window, 'Position');
height = pos(4);

%Set y-location for Axis_Left, Axis_Center and Axis_Right (pixels)
y = [138; 129; 138];

%Fail-safe: return when height falls below y
if height<max(y), return; end

%Modify positions to comply
y_orig = positions(:,2);
positions(:,2) = y/height;
positions(:,4) = positions(:,4)+y_orig-positions(:,2);

%Set axes size
set(ax(1), 'Position', positions(1,:));
set(ax(2), 'OuterPosition', positions(2,:));
set(ax(3), 'Position', positions(3,:));



%-------------------%
% Utility functions %
%-------------------%


%Utility function, returns visible axes
function ax = VisibleAxes(handles)

%Get axes handles
ax = [handles.Left_Axis; handles.Center_Axis; handles.Right_Axis];

%Limit ax to axes currently visible
ax = ax(strcmp(get(ax, 'Visible'), 'on'));



%Utility function, toggles control panel visibility.
function TogglePanel(handles)

%Get handles controls and panels
h = [handles.Walk; get(handles.Walk, 'Children'); ...
     handles.Orientation; get(handles.Orientation, 'Children'); ...
     handles.Zoom; get(handles.Zoom, 'Children'); ...
     handles.Stereo; get(handles.Stereo, 'Children'); ...
     handles.Perspective; handles.Text; handles.IDC; handles.Backpanel];

%Get handles of panels only
h_p = [handles.Walk; handles.Orientation; handles.Zoom; handles.Stereo];

%Toggle visibility
if any(strcmp(get(h, 'Visible'), 'off'))
    set(h, 'Visible', 'on');
    set(h_p, 'BackgroundColor', 'default', 'ForegroundColor', 'default');
else
    set(h, 'Visible', 'off');
    color = get(handles.Window, 'Color');
    set(h_p, 'BackgroundColor', color, 'ForegroundColor', color);
end



%Utility function, separates the cameras
function Stereo(ax_lefteye, ax_righteye, idc)

%Get normalised camera position, camera up vector and target
aspect = get(ax_lefteye, 'DataAspectRatio');
pos = campos(ax_lefteye)./aspect;
tar = camtarget(ax_lefteye)./aspect;
upv = camup(ax_lefteye)./aspect;
d = sqrt(sum((pos-tar).^2));

%Calculate the position of the two cameras <unitv>
l = d/idc;
pos1 = pos+l/2*unitv(cross(pos-tar, upv));
pos2 = pos-l/2*unitv(cross(pos-tar, upv));

%Make sure camera doesn't "walk" after multiple switches
pos1 = tar+d*unitv(pos1-tar);
pos2 = tar+d*unitv(pos2-tar);

%Set up new positions
campos(ax_lefteye, pos1.*aspect);
campos(ax_righteye, pos2.*aspect);
camtarget(ax_lefteye, tar.*aspect);
camtarget(ax_righteye, tar.*aspect);

%Reset camera up vector
set([ax_lefteye; ax_righteye], 'CameraUpVectorMode', 'auto');



%Utility function, normalises the input
function u = unitv(v)

%Divide v by its norm
u = v/(sqrt(sum(v.^2)));



% %Original GUIDE help text (wb20060310), begin:
% 
% GUI_MESHEXPLORER M-file for GUI_MeshExplorer.fig
%      GUI_MESHEXPLORER, by itself, creates a new GUI_MESHEXPLORER.
%
%      H = GUI_MESHEXPLORER returns the handle to a new GUI_MESHEXPLORER or the handle to
%      the existing singleton*.
%
%      GUI_MESHEXPLORER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_MESHEXPLORER.M with the given input arguments.
%
%      GUI_MESHEXPLORER('Property','Value',...) creates a new GUI_MESHEXPLORER or raises the
%      existing singleton*.  Starting from the left_axis, property value pairs are
%      applied to the GUI before GUI_MeshExplorer_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_MeshExplorer_OpeningFcn via varargin.
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% 
% %Original GUIDE help text (wb20060310), end.