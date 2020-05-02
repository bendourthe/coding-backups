function varargout = GUI_MainInterface(varargin)
% GUI_MAININTERFACE M-file for GUI_MainInterface.fig
%      GUI_MAININTERFACE, by itself, creates a new GUI_MAININTERFACE or raises the existing
%      singleton*.
%
%      H = GUI_MAININTERFACE returns the handle to a new GUI_MAININTERFACE or the handle to
%      the existing singleton*.
%
%      GUI_MAININTERFACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_MAININTERFACE.M with the given input arguments.
%
%      GUI_MAININTERFACE('Property','Value',...) creates a new GUI_MAININTERFACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_MainInterface_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_MainInterface_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_MainInterface

% Last Modified by GUIDE v2.5 18-Aug-2006 10:54:13



%----------------%
% Initialisation %
%----------------%


%Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_MainInterface_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_MainInterface_OutputFcn, ...
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



%Executed just before GUI_MainInterface is made visible
function GUI_MainInterface_OpeningFcn(src, evt, h, varargin)

%Update handles structure
guidata(src, h);



%Outputs from this function are returned to the command line
function varargout = GUI_MainInterface_OutputFcn(src, evt, h) 

%Get default command line output from handles structure
varargout{1} = h.MainFigure;



%-----------%
% Callbacks %
%-----------%


%Callback, executed when MainFigure is resized
function MainFigure_ResizeFcn(src, evt, h)

%Set parameter
border = 10; %Border between panels

%Make sure units are set to pixels
set([h.MainFigure; h.TextPanel; h.AxesPanel; h.ControlPanel], 'Units', 'pixels');

%Store positions of objects
pos_tp = get(h.TextPanel, 'Position');
pos_ap = get(h.AxesPanel, 'Position');
pos_cp = get(h.ControlPanel, 'Position');
pos_mf = get(h.MainFigure, 'Position');

%Keep TextPanel on the bottom left, adapt width
pos_tp(1) = border+1;
pos_tp(2) = border+1;
pos_tp(3) = max(pos_mf(3)-pos_cp(3)-3*border, 1);
set(h.TextPanel, 'Position', pos_tp);
pos_tp = get(h.TextPanel, 'Position'); %In case TextPanel's ResizeFcn modified its position

%Keep AxesPanel on the left and on top of TextPanel, adapt width and height
pos_ap(1) = border+1;
pos_ap(2) = pos_tp(2)+pos_tp(4)+border;
pos_ap(3) = pos_tp(3);
pos_ap(4) = max(pos_mf(4)-pos_ap(2)+1-border, 1);
set(h.AxesPanel, 'Position', pos_ap);
pos_ap = get(h.AxesPanel, 'Position'); %In case AxesPanel's ResizeFcn modified its position

%Make sure TextPanel's width follows that of AxesPanel
if pos_ap(3)>pos_tp(3)
    pos_tp(3) = pos_ap(3);
    set(h.TextPanel, 'Position', pos_tp);
end

%Keep ControlPanel on the bottom right, adapt height
pos_cp(1) = pos_tp(1)+pos_tp(3)+border;
pos_cp(2) = border+1;
pos_cp(4) = pos_tp(4)+border+pos_ap(4);
set(h.ControlPanel, 'Position', pos_cp);
pos_cp = get(h.ControlPanel, 'Position');

%Make sure AxesPanel's top edge follows that of ControlPanel
if pos_cp(2)+pos_cp(4)>pos_ap(2)+pos_ap(4)
    pos_ap(4) = pos_cp(2)+pos_cp(4)-pos_ap(2);
    set(h.AxesPanel, 'Position', pos_ap);
end

%Adapt MainFigure's size to panels if necessary
pos_mf_min(1) = pos_cp(1)-1+pos_cp(3)+border;
pos_mf_min(2) = pos_ap(2)-1+pos_ap(4)+border;
if any(pos_mf_min>pos_mf(3:4))
    pos_mf(3:4) = pos_mf_min;
    set(h.MainFigure, 'Position', pos_mf);
end



%Callback, executed when TextPanel is resized
function TextPanel_ResizeFcn(src, evt, h)

%Set parameters
border = 2; %Border between Text and TextPanel
w_min = 50; %Minimum width of Text
h_min = 10; %Minimum height of Text

%Make sure units are set to pixels
set([0; h.MainFigure; h.TextPanel; h.Text], 'Units', 'pixels');

%Store positions of objects
pos_tp = get(h.TextPanel, 'Position');
pos_tx = get(h.Text, 'Position');

%Determine difference between coordinates of Text and Textpanel
offset = hgconvertunits(h.MainFigure, [0 0 1 1], 'normalized', 'pixels', h.TextPanel);
offset = offset-[1 1 pos_tp(3:4)];

%Resize Text
pos_tx(1) = offset(1)+border+1;
pos_tx(2) = offset(2)+border+1;
pos_tx(3) = max(offset(3)-offset(1)+pos_tp(3)-2*border, w_min);
pos_tx(4) = max(offset(4)-offset(2)+pos_tp(4)-2*border, h_min);
set(h.Text, 'Position', pos_tx);

%Adapt TextPanel's size to Text if necessary
pos_tp_min = offset(1:2)+pos_tx(3:4)+2*border-offset(3:4);
if any(pos_tp_min>pos_tp(3:4))
    pos_tp(3:4) = pos_tp_min;
    set(h.TextPanel, 'Position', pos_tp);
end



%Callback, executed when the user clicks on Text
function Text_ButtonDownFcn(src, evt, h)

%Set parameter
fraction = 0.5; %Fraction of figure height used when enlarging

%Store height of figure space
space = get(h.MainFigure, 'Position');
space = space(4);

%Determine size of TextPanel
pos = get(h.TextPanel, 'Position');

%If application data is present, shrink; otherwise, enlarge
if isappdata(h.TextPanel, mfilename)
    
    %Set height of TextPanel to stored data
    pos(4) = getappdata(h.TextPanel, mfilename);
    
    %Remove application data
    rmappdata(h.TextPanel, mfilename);
    
else
    
    %Fail-safe: return if figure height is too small
    if fraction*space<=pos(2)+pos(4)-1
        return
    end
    
    %Store original TextPanel height in application data
    setappdata(h.TextPanel, mfilename, pos(4));
    
    %Modify height of TextPanel
    pos(4) = fix(fraction*space-pos(2)+1);
end

%Update TextPanel position <MainFigure_ResizeFcn>
set(h.TextPanel, 'Position', pos);
MainFigure_ResizeFcn(src, evt, h);