function [selected, ax] = GUI_SelectShells(varargin)
%GUI_SelectShells  Select shells using GUI.
%wb20060620
%
%   Syntax:
%    [selected, ax] = ...
%        GUI_SelectShells(F, V, colors, colors_sel, min_sel, max_sel)
%    [selected, ax] = ...
%        GUI_SelectShells(h, F, V, colors, colors_sel, min_sel, max_sel)
%
%   Input:
%    h:          Handle to the figure or axes on which the mesh is to be
%                displayed. If h is an axes handle, the mesh will be added
%                to this axes. If h is a figure handle, the mesh will be
%                displayed in a new axes on h. If h is 0, a new figure will
%                be created.
%    F:          N-by-1 cell array containing separate F-matrices for each
%                shell. Each F-matrix is an N-by-3 array containing indices
%                into the corresponding V. Each row represents a triangle,
%                each element is a link to a vertex in V.
%    V:          N-by-1 cell array containing separate V-matrices for each
%                shell. Each V-matrix is an N-by-3 array containing vertex
%                coordinates. Each row represents a vertex; the first,
%                second and third columns represent X-, Y- and Z-
%                coordinates respectively.
%    colors:     N-by-3 array containing color definitions. Each row
%                corresponds with a shell in F and V. The first, second and
%                third columns represent red, green and blue values
%                respectively. Each element lies between 0 and 1. If set to
%                an empty matrix, colors will be assigned automatically.
%                Optional, defaults to an empty matrix.
%    colors_sel: N-by-3 array containing color definitions to be used for
%                selected shells. Each row corresponds with a shell in F
%                and V. The first, second and third columns represent red,
%                green and blue values respectively. Each element lies
%                between 0 and 1. If set to an empty matrix, colors will be
%                assigned automatically. Optional, defaults to an empty
%                matrix.
%    min_sel:    The minimum number of selected shells. The user cannot
%                close the window until this number of selected shells has
%                been reached. Optional, defaults to 0.
%    max_sel:    The maximum number of selected shells. The user cannot
%                select any additional shells when this number of selected
%                shells has been reached. Optional, defaults to Inf.
%
%   Output:
%    selected: Column vector of logicals indicating whether or not the
%              corresponding shell was selected.
%    ax:       Handle to the axes on which the mesh is displayed.
%
%   Effect: This function will display the mesh defined by F and V on a new
%   figure or on the object represented by h, depending on which syntax is
%   used. The user can select one or more of the shells in the mesh. The
%   selection ends when the enter key is pressed, provided that the minimum
%   number of selected shells, specified by min_sel, has been reached. If
%   the axes handle is not requested by the calling function, the window
%   will close when the user ends the selection; otherwise, the window will
%   stay open.
%
%   Dependencies: GUI_DistributeColors.m
%                 GUI_PlotShells.m
%                 GUI_VisualisationKeyPress.m
%
%   Known parents: TRI_MeshCutter.m
%                  TRI_RegionCutter.m

%Created on 19/04/2006 by Ward Bartels.
%WB, 03/05/2006: Figure handles are now accepted as input arguments.
%WB, 20/06/2006: Closing the window is now possible (default output
%                arguments will be passed).
%Stabile, fully functional.



%----------------%
% Initialisation %
%----------------%


%Handle input
if numel(varargin{1})==1 && ishandle(varargin{1})
    [h, F, V] = varargin{1:3}; %h may be 0!
    argin = varargin(4:end);
    nin = nargin-3;
else
    [F, V] = varargin{1:2};
    h = 0;
    argin = varargin(3:end);
    nin = nargin-2;
end

%Set defaults for colors, colors_sel, min_sel and max_sel
defaults = {[] [] 0 Inf};
argin(4:-1:nin+1) = defaults(4:-1:nin+1);

%Transfer inputs
[colors, colors_sel, min_sel, max_sel] = argin{:};

%Fail-safe, prevents the user from getting stuck in this window
if min_sel>max_sel
    error([mfilename ': The minimum number of selected shells cannot be greater than the maximum number of selected shells.']);
end
if min_sel>length(F)
    error([mfilename ': The minimum number of selected shells cannot be greater than the number of shells provided.']);
end

%Automatically assign colors if necessary <<GUI_DistributeColors.m>>
if isempty(colors)
    colors = GUI_DistributeColors(length(F));
elseif size(colors, 1)~=length(F)
    colors = colors(ones(length(F), 1),:);
end

%Automatically assign colors_sel if necessary
if isempty(colors_sel)
    colors_sel = 1-colors; %Invert colors
elseif size(colors_sel, 1)~=length(F)
    colors_sel = colors_sel(ones(length(F), 1),:); %Replicate single color
end

%Plot the mesh and enable selection <SelectShell> <<GUI_PlotShells.m>>
[obj, li, ax] = GUI_PlotShells(h, F, V, colors, 0, [], false);

%Turn off figure interactive properties
fig = ancestor(ax, 'figure');
state = uiclearmode(fig, 'docontext');

%Save variables to application data
ad = struct('ax', ax, 'li', li, 'obj', obj, 'selected', [], 'colors', colors, ...
            'colors_sel', colors_sel, 'min_sel', min_sel, 'max_sel', max_sel);
setappdata(fig, mfilename, ad);

%Set figure and object callbacks <SelectShell> <KeyPress>
%                                <<GUI_VisualisationKeyPress.m>>
set(obj, 'ButtonDownFcn', @SelectShell);
set(fig, 'KeyPressFcn', {@GUI_VisualisationKeyPress, ax, @KeyPress});

%Put figure on foreground
figure(fig);

%Suspend execution
uiwait(fig);

%If user closed the window, return
if ~ishandle(fig)
    selected = false(0, 0);
    ax = NaN;
    return
end

%Get application data
ad = getappdata(fig, mfilename);
rmappdata(fig, mfilename);

%Reset object properties and disable selection
set(ad.obj, 'ButtonDownFcn', '');
set(ad.obj, {'FaceColor'}, num2cell(ad.colors, 2));
uiclearmode(fig, 'docontext');
uirestore(state);

%Determine which objects were selected
selected = ismember(ad.obj, ad.selected);

%Close the window if axes handle isn't requested in the output
if nargout<2, close(fig); end



%-----------%
% Callbacks %
%-----------%


%Callback, executed when the user clicks on a shell
function SelectShell(src, evt)

%Get application data
fig = gcbf;
ad = getappdata(fig, mfilename);

%Determine whether or not the selected object is already selected
ind = ad.selected==src;
if any(ind)                    %Selected
    
    %Reset color of selected shell
    set(src, 'FaceColor', ad.colors(ad.obj==src,:));
    
    %Unmark hgtransform object as selected
    ad.selected(ind,:) = [];
    
elseif size(ind, 1)<ad.max_sel %Not yet selected and maximum not reached
    
    %Set color of selected shell
    set(src, 'FaceColor', ad.colors_sel(ad.obj==src,:));
    
    %Mark object as selected
    ad.selected(end+1,1) = src;
end

%Update application data
setappdata(fig, mfilename, ad);



%Callback, executed when a key is pressed in the figure
function KeyPress(fig, evt)

%Get application data
ad = getappdata(fig, mfilename);

%If enter was pressed, allow main function to continue
if strcmp(evt.Key, 'return')
    
    %If enough shells selected, let main function proceed; otherwise, warn
    if size(ad.selected, 1)>=ad.min_sel
        uiresume(fig);
    else
        warndlg(['Please select at least ' num2str(ad.min_sel) ' shell(s).'], ...
                 'Not enough shells selected', 'modal');
        beep;
    end
end