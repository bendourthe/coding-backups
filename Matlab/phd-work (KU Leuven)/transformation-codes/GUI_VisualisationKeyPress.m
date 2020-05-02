function GUI_VisualisationKeyPress(src, evt, ax, varargin)
%GUI_VisualisationKeyPress  Standard key press callback for visualisations.
%
%   Syntax:
%    GUI_VisualisationKeyPress(src, evt, ax, args{:})
%
%   Input:
%    src:  Handle to the object generating the callback. Automatically
%          passed by Matlab when function handle syntax is used.
%    evt:  Event data structure. Automatically passed by Matlab when
%          function handle syntax is used.
%    ax:   Column vector of handles to the axes on which this function will
%          take effect. If ax is empty, the axes under src will be used.
%          Optional, defaults to an empty matrix.
%    args: Contents of a callback property that will be activated if this
%          function does not recognise the event.
%
%   Effect: This function is a standard key press callback intended to
%   provide additional functionality to visualisation windows. When this
%   function is set as a figure's KeyPressFcn, the following key
%   combinations may be used:
%      - Space:       Toggle last used rotation in Visualisation toolbar
%      - Up arrow:    Tilt mesh upward
%      - Down arrow:  Tilt mesh downward
%      - Left arrow:  Rotate mesh towards the left
%      - Right arrow: Rotate mesh towards the right
%      - Page Up:     Zoom in
%      - Page Down:   Zoom out
%      - Home:        Reset roll and zoom
%      - L:           Line up lights with the camera
%      - S:           View current situation with Mesh Explorer
%   If any other keys are used, the event is passed on to the callback
%   defined by args.
%
%   Example:
%    set(gcf, 'KeyPressFcn', {@GUI_VisualisationKeyPress, gca});
%
%   Dependencies: GUI_MeshExplorer.m
%                 GUI_RunCallback.m
%
%   Known parents: GUI_VisualisationUI.m
%                  GUI_IndicateMultiPlane.m
%                  GUI_SimulateSurgery.m
%                  GUI_SelectShells.m
%                  TRI_RegionCutter.m

%Created on 19/04/2006 by Ward Bartels.
%WB, 03/05/2006: Deleting the axes will no longer cause errors; axes handle
%                can now be determined automatically.
%WB, 24/05/2006: Added zoom and roll functionality.
%WB, 05/07/2006: Added axes linking.
%WB, 25/09/2006: GUI_RunCallback is now used to pass on the event.
%Stabile, fully functional.


%Check if modifiers were used
passon = ~isempty(evt.Modifier);

%Get axes handle if necessary
if ~passon && (nargin<3 || isempty(ax))
    ax = findobj(src, 'Type', 'axes');
elseif ~passon
    ax = ax(ishandle(ax));
end

%Pass on event if ax is gone
passon = passon || isempty(ax);

%Determine which key was pressed <<GUI_MeshExplorer.m>>
if ~passon
    switch evt.Key
        case 'space'      %Space - activate last rotation mode
            tb = findobj(src, 'Tag', 'VisToolbar');
            if isempty(tb), return; end
            button = getappdata(tb(1), 'LastViewButton');
            if isempty(button), return; end
            if strcmp(get(button, 'State'), 'on')
                set(button, 'State', 'off');
            else
                set(button, 'State', 'on');
            end
        case 'home'       %Home - reset view
            set(ax, 'CameraTargetMode', 'auto');
            set(ax, 'CameraPositionMode', 'auto');
            set(ax, 'CameraViewAngleMode', 'auto');
            set(ax, 'CameraUpVectorMode', 'auto');
        case 'uparrow'    %Up arrow - decrease elevation
            arrayfun(@(x) camorbit(x, 0, -5), ax);
        case 'downarrow'  %Down arrow - increase elevation
            arrayfun(@(x) camorbit(x, 0, 5), ax);
        case 'rightarrow' %Right arrow - decrease azimuth
            arrayfun(@(x) camorbit(x, -5, 0), ax);
        case 'leftarrow'  %Left arrow - increase azimuth
            arrayfun(@(x) camorbit(x, 5, 0), ax);
        case 'end'        %Numpad 1 - roll camera clockwise
            arrayfun(@(x) camroll(x, 5), ax);
        case 'insert'     %Numpad 0 - roll camera counterclockwise
            arrayfun(@(x) camroll(x, -5), ax);
        case 'pageup'     %Page Up - zoom in
            arrayfun(@(x) camzoom(x, 1.1), ax);
        case 'pagedown'   %Page Down - zoom out
            arrayfun(@(x) camzoom(x, 0.9), ax);
        case 's',         %S - start Mesh Explorer
            uiwait(GUI_MeshExplorer(gca));
        case 'l',         %L - move lights
            for ind = 1:numel(ax)
                li = get(ax(ind), 'Children');
                li = li(strcmp(get(li, 'Type'), 'light'));
                if length(li)>=1
                    camlight(li(1), 'headlight', 'infinite');
                end
                if length(li)>=2
                    camlight(li(2), 0, 180, 'infinite');
                end
            end
        otherwise,        %Unknown - pass event to fun
            passon = true;
    end
end

%Pass event to fun <<GUI_RunCallback.m>>
if passon && nargin>=4
    GUI_RunCallback(src, evt, varargin{:});
end