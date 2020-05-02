function varargout = mech_import_stl(varargin)
% Copyright 2011 The MathWorks, Inc.

  % MECH_IMPORT_STL MATLAB code for mech_import_stl.fig
  %      MECH_IMPORT_STL, by itself, creates a new MECH_IMPORT_STL or raises the existing
  %      singleton*.
  %
  %      H = MECH_IMPORT_STL returns the handle to a new MECH_IMPORT_STL or the handle to
  %      the existing singleton*.
  %
  %      MECH_IMPORT_STL('CALLBACK',hObject,eventData,handles,...) calls the local
  %      function named CALLBACK in MECH_IMPORT_STL.M with the given input arguments.
  %
  %      MECH_IMPORT_STL('Property','Value',...) creates a new MECH_IMPORT_STL or raises the
  %      existing singleton*.  Starting from the left, property value pairs are
  %      applied to the GUI before mech_import_stl_OpeningFcn gets called.  An
  %      unrecognized property name or invalid value makes property application
  %      stop.  All inputs are passed to mech_import_stl_OpeningFcn via varargin.
  %
  %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
  %      instance to run (singleton)".
  %
  % See also: GUIDE, GUIDATA, GUIHANDLES
  
  % Edit the above text to modify the response to help mech_import_stl
  
  % Last Modified by GUIDE v2.5 25-Mar-2011 08:59:33
  
  % Begin initialization code - DO NOT EDIT
  gui_Singleton = 1;
  gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @mech_import_stl_OpeningFcn, ...
    'gui_OutputFcn',  @mech_import_stl_OutputFcn, ...
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
  % End initialization code - DO NOT EDIT
  
  
  % --- Executes just before mech_import_stl is made visible.
function mech_import_stl_OpeningFcn(hObject, eventdata, handles, varargin)
  % This function has no output args, see OutputFcn.
  % hObject    handle to figure
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  % varargin   command line arguments to mech_import_stl (see VARARGIN)
  
  % Choose default command line output for mech_import_stl
  handles.output = hObject;
  
  % Update handles structure
  guidata(hObject, handles);
  
  % UIWAIT makes mech_import_stl wait for user response (see UIRESUME)
  % uiwait(handles.figure1);
  
  
  % --- Outputs from this function are returned to the command line.
function varargout = mech_import_stl_OutputFcn(hObject, eventdata, handles)
  % varargout  cell array for returning output args (see VARARGOUT);
  % hObject    handle to figure
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  
  % Get default command line output from handles structure
  varargout{1} = handles.output;
  
  
  % --- Executes during object creation, after setting all properties.
function file_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
  % hObject    handle to file (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    empty - handles not created until after all CreateFcns called
  
  % Hint: edit controls usually have a white background on Windows.
  %       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end
  
  
  % --- Executes on button press in pb_browse_stl.
function pb_browse_stl_Callback(hObject, eventdata, handles)
  % hObject    handle to pb_browse_stl (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  filename = stl2m();
  handles.filename = filename;
  set(handles.file,'String',filename(1:end-2));
  set(handles.pb_create,'Enable','on');
  % Update handles structure
  guidata(hObject, handles);
  
  
  % --- Executes on button press in pb_create.
function pb_create_Callback(hObject, eventdata, handles)
  % hObject    handle to pb_create (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  contents = cellstr(get(handles.pu_material,'String'));
  material = contents{get(handles.pu_material,'Value')};
  density = handles.d(material);
  
  [M,Xg,J] = stl2sm(handles.filename(1:end-2),density);
  
  load_system('mblibv1');
  sys='untitled';
  try
    open_system(sys);
  catch me %#ok
    new_system(sys);
    open_system(sys);
  end
  try
    bh=add_block('mblibv1/Bodies/Body',[sys '/' handles.filename(1:end-2)]);
    
    set(bh,'Mass',num2str(M));
    set(bh,'Inertia',mat2str(J));
    s=get(bh,'CG');
    s=strrep(s,'WORLD','CS1');
    s=strrep(s,'CG$[0 0 0]',['CG$' mat2str(Xg')]);
    set(bh,'CG',s);
    s=get(bh,'WorkingFrames');
    s=strrep(s,'CS1$[0 0 0]$CG$CG','CS1$[0 0 0]$Adjoining$Adjoining');
    set(bh,'WorkingFrames',s);
     
    %Visuallization
    set(bh,'AttachedToCS','CS1');
    set(bh,'GraphicsMode','GFXFILE');
    set(bh,'GraphicsFileName',[handles.filename(1:end-2) '.stl']);
    set_param(sys,'VisOnUpdateDiagram','on');
    set_param(sys,'VisDuringSimulation','on');
    
    %Apply changes
    mech_support('LoadFcn',bh);
  catch me
    me %#ok
  end
  
  
  % --- Executes during object creation, after setting all properties.
function pu_material_CreateFcn(hObject, eventdata, handles)
  % hObject    handle to pu_material (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    empty - handles not created until after all CreateFcns called
  
  % Hint: popupmenu controls usually have a white background on Windows.
  %       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end
  
  % Init list box
  d=density();
  set(hObject,'String',d.keys);
  handles.d = d;

  % Update handles structure
  guidata(hObject, handles);
  


% --- Executes on selection change in pu_material.
function pu_material_Callback(hObject, eventdata, handles)
% hObject    handle to pu_material (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pu_material contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pu_material
