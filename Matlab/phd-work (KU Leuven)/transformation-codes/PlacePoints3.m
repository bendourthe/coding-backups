function [ad] =  PlacePoints3(F, V, F_surf, V_surf, title)

%PlacePoints  Interactively cuts multiple regions from a mesh based on landmarks.
%
%   Syntax:
%    [F_part, V_part, ad] = ...
%        TRI_MarkerRegionCutter(F, V, markers, select, colors, ad_load)
%
%   Input:
%    F:      N-by-3 array containing indices into V. Each row represents a
%            triangle, each element is a link to a vertex in V. If the mesh
%            contains multiple shells, the F-matrices should be placed in
%            an N-by-1 cell array. 
%    V:      N-by-3 array containing vertex coordinates. Each row
%            represents a vertex; the first, second and third columns
%            represent X-, Y-� and Z-coordinates respectively. If the mesh
%            contains multiple shells, the V-matrices should be placed in
%            an N-by-1 cell array. If it is empty, V will be obtained from
%            an STL-file selected through a GUI. Optional, defaults to an
%            empty cell array.
%    colors: N-by-3 array containing color definitions. Each row
%            corresponds with a shell in F and V. The first, second and
%            third columns represent red, green and blue values
%            respectively. Each element lies between 0 and 1. If set to an
%            empty matrix, colors will be assigned automatically. Optional,
%            defaults to an empty matrix.
%   
%
%   Output:
%    F:     N-by-3 array containing indices into V. Each row represents a
%           triangle in the shell cut out by the user, each element is a
%           link to a vertex in V.
%    V:     N-by-3 array containing vertex coordinates. Each row represents
%           a vertex; the first, second and third columns represent X-, Y-
%           and Z-coordinates respectively.
%    F_rem: N-by-3 array containing indices into V. Each row represents a
%           triangle in the remainder of the shell cut out by the user,
%           each element is a link to a vertex in V.
%    V_rem: N-by-3 array containing vertex coordinates. Each row represents
%           a vertex; the first, second and third columns represent X-, Y-
%           and Z-coordinates respectively.
%    shell: Row index in the original F and V of the shell that was cut.
%    ax:    Handle to the axes on which the mesh is displayed.
%

%   Dependencies:   CreateColors.m
%                   GUI_VisualisationUI.m
%                   GUI_PlotShells
%                   GUI_VisualisationKeyPress

% FK, 18/4/2013: added dependencies
%----------------%
% Initialisation %
%----------------%

nofig = 1;
F_m = F;
V_m = V;
%Plot the mesh <<GUI_PlotShells.m>>
ax = zeros(nofig,1);
obj = cell(nofig,1);
for i = 1: nofig
    [obj_tmp, ign, ax(i)] = GUI_PlotShells(F{i}, V{i}, CreateColors(size(F{i},1)));
    obj(i) = {obj_tmp};
    GUI_VisualisationUI(ancestor(ax(i), 'figure'), true, ax(i), true, true);
    if nargin==5
    set(gcf,'Name',title)
    end
end

% %Plot invisible meshes on the second image
obj_inv = cell(nofig,1);
% for i = setdiff(1:nofig,2)
%     obj_tmp = GUI_PlotShells(ax(2),F{i}, V{i}, CreateColors(size(F{i},1)), 0, [], true);
%     obj_inv(i) = {obj_tmp};
%     set(obj_tmp,'Visible','off');
% end

%Turn off figure interactive properties
fig = zeros(nofig,1);
state = cell(nofig,1);
for i = 1:nofig
    fig(i) = ancestor(ax(i), 'figure');
    state(i) = {uiclearmode(fig(i), 'docontext')};
end

%Set hold state
hstate = false(nofig,1);
for i = 1:nofig
    hstate(i) = ishold(ax(i));
    hold(ax(i), 'on');
end

%Put first figure on foreground
figure(fig(1));

%Get a length scale
scale = zeros(nofig,1);
for i = 1:nofig
    scale(i) = diff(get(ax(i), 'XLim'));
end

% Parts
T =  cell(nofig,1);
for i = 1:nofig
    T(i) = {repmat(eye(4),[1,1,size(obj{i},1)])};
end


%Initialise variables and save to application data to the first figure
ad = struct('ax', ax, 'scale', scale, 'pl', [], 'F_surf', {F_surf}, 'V_surf', {V_surf}, ...
            'pts', {cell(0, 1)}, 'selected', 0, 'figs', fig, 'T', {T}, 'V_surf_old', {V_surf}, ...
            'F_m', {F_m}, 'V_m', {V_m}, 'obj_inv', {obj_inv}, 'shell', {cell(0,1)}, ...
            'F', {F}, 'V', {V}, 'points', {cell(0,1)}, 'obj', {obj}, ...
            'F_orig', {F}, 'V_orig', {V}, 'F_m_orig', {F_m}, 'V_m_orig', {V_m}, ...
            'copy', false, 'move', false, 'cur_fig', 1, 'nofig', nofig);
setappdata(fig(1), mfilename, ad);

%Set callbacks <ButtonPress> <KeyPress> <<GUI_VisualisationKeyPress.m>>
set(ad.figs(ad.cur_fig), 'WindowButtonDownFcn', {@ButtonPress fig}, ...
         'KeyPressFcn', {@GUI_VisualisationKeyPress ax(i) @KeyPress fig});

%Suspend execution
uiwait(fig(1));

%If user closed the window, return
for i = 1:nofig
    if ~ishandle(fig(i))
        [pts] = deal(zeros(0, 3)); %???
        ad = deal(0); %???
        return
    end
end

%Get application data
ad = getappdata(fig(1), mfilename);

% set(ad.pts{end,1},'Visible','off');

%Reset figure properties
for i = 1: nofig
    uiclearmode(fig(i), 'docontext');
    uirestore(state{i});
    if ~hstate(i), hold(ax(i), 'off');
end
end



%-----------%
% Callbacks %
%-----------%


%Callback, executed when a mouse button is pressed in the figure
function ButtonPress(fig, evt, figs)

%Get application data
ad = getappdata(figs(1), mfilename);

%Get selected point
selection = get(ad.ax(ad.cur_fig(end)), 'CurrentPoint');

% if size(ad.pts,1)
%     disp('Before AddPoint');
%     disp('pts');
%     disp(ad.pts{end,1});
%     disp('points');
%     disp(ad.points{end,1});
%     disp('selected');
%     disp(ad.selected(end));
%     disp('shell');
%     disp(ad.shell{end,1});
% end

%Close the cut if the right button was pressed
if strcmp(get(fig, 'SelectionType'), 'alt')
    %change figure
    ad = ChangeFigure(ad);
elseif ad.move
    ad = MovePoint(ad,selection);
else
    %Cut the mesh <AddCut>
    ad = AddPoint(ad, selection, false);
end
ad = RefreshSelection(ad);
% 
% if size(ad.pts,1)
%     disp('After AddPoint');
%     disp('pts');
%     disp(ad.pts{end,1});
%     disp('points');
%     disp(ad.points{end,1});
%     disp('selected');
%     disp(ad.selected(end));
%     disp('shell');
%     disp(ad.shell{end,1});
% end

%Store application data
setappdata(figs(1), mfilename, ad);

%Callback, executed when a key is pressed in the figure
function KeyPress(fig, evt, figs)

%Get application data
ad = getappdata(figs(1), mfilename);

%Determine which key was pressed <DeleteLastPoly>
switch evt.Key
    case 'return' %Enter - end selection (if path is closed)
        uiresume(figs(1));
            return
%             else
%                 warndlg('Please select a closed cutting path.', ...
%                         'Path not closed', 'modal');
%                 beep;
%             end
    case 'delete' %Delete - delete selected point
        ad = DeletePoint(ad);
        ad = RefreshSelection(ad);
    case 't' %Delete - delete selected point
        ad = DeletePoint(ad);
        ad = RefreshSelection(ad);
    case 'r' %Register
        if strcmp(evt.Modifier,'control')
            ad = Registration(ad,'parts ICP');
        elseif strcmp(evt.Modifier,'alt')
            ad = Registration(ad,'parts paired');
        else
            ad = Registration(ad,'full');
        end
        ad = RefreshSelection(ad);
    case 'p'      %cycle through contourselection forwards
        selected = ad.selected(end);
        pts = ad.pts{end,1};
        if selected < size(pts,1)
            ad.selected(end) = selected +1;  
        else
            ad.selected(end) = 1;
        end
            ad = RefreshSelection(ad);
    case 'm'      %cycle through contourselection backwards
        selected = ad.selected(end);
        pts = ad.pts{end,1};
        if selected > 1
            ad.selected(end) = selected - 1;
        else
            ad.selected(end) = size(pts,1);
        end
            ad = RefreshSelection(ad);
    case 'o' % toggle overlay
        ad = ToggleOverlay(ad);
    case 'i' % toggle facealpha
        ad = ToggleFaceAlpha(ad);
    case 'z'
        if strcmp(evt.Modifier,'control')
            ad = UndoAction(ad);
        end
    case 'n' % move a point
        ad.move = ~ad.move;
        ad.copy = false;
end

%Update application data
setappdata(figs(1), mfilename, ad);

%-------------------%
% Utility functions %
%-------------------%

%Utility function, adds a point
function ad = AddPoint(ad, selection, auto)
if size(ad.points,1) == 0;          % first data
    points = zeros(0,3*ad.nofig);
    pts = zeros(0,1);
    selected = 0;
    shell = zeros(0,1);
else                                % next points
    pts = ad.pts{end,:};
    shell = ad.shell{end,:};
    points = ad.points{end,:};
    selected = ad.selected(end);
end
cur_fig = ad.cur_fig(end);
% Find the position of the new point
% direction = diff(selection,1,1);
F_tmp = ad.F_m_orig{cur_fig};
V_tmp = ad.V_m_orig{cur_fig};
% point on mesh
[point, triangles] = TRI_IntersectWithLine(F_tmp, V_tmp, selection);
% find point closest to the camera
if auto
    [ignoble, ind] = min(sum((point-ones(size(point, 1), 1)*selection(1,:)).^2, 2));
else
    [ignoble, ind] = min(sum((point-ones(size(point, 1), 1)*campos(ad.ax(cur_fig))).^2, 2));
end
if ~isempty(ind) % closest point on the intersected triangle's edge
    point = point(ind,:);
    point = ClosestEdgePoint(point,V_tmp(F_tmp(triangles(ind),:),:));
else % no intersect with the mesh
    return
end
pts_new = plot3(point(1), point(2), point(3), 'k.', 'Parent', ad.ax(cur_fig)).';

% Check which shell the point belongs to
V_tmp2 = ad.V_orig{cur_fig};
shell_new = 0;
if iscell(V_tmp2)
    for i = 1: size(V_tmp2,1);
        V_tmp3 = V_tmp2{i};
        if ismember(V_tmp(F_tmp(triangles(ind),1),:), V_tmp3 ,'rows') %any(all(V_tmp3-ones(size(V_tmp3,1),1)*point,2))
            shell_new = i;
        end
    end
end

% Check where to add the point and add them
if ad.move && selected ~= 0
    %replace point
    points(selected,3*(cur_fig-1)+1:3*cur_fig) = point;
    pts(selected,ad.cur_fig) = pts_new;
elseif ~size(pts,1)
    %new point
    points = [points; NaN(1,3*ad.nofig)];
    points(selected+1,3*(cur_fig-1)+1:3*cur_fig) = point;
    pts = [pts; NaN(1,ad.nofig)];
    pts(selected+1,cur_fig) = pts_new;
    shell = [shell; NaN(1,ad.nofig)];
    shell(selected+1,cur_fig) = shell_new;
    selected = selected+1;
elseif selected <= size(pts,1)
    if isnan(pts(selected,cur_fig))
        points(selected,3*(cur_fig-1)+1:3*cur_fig) = point;
        pts(selected,cur_fig) = pts_new;
        shell(selected,cur_fig) = shell_new;
    elseif selected < size(pts,1)
        %intermediate point
        points = [points(1:selected,:); NaN(1,3*ad.nofig); points(selected+1:size(points,1),:)];
        points(selected+1,3*(cur_fig-1)+1:3*cur_fig) = point;
        pts = [pts(1:selected,:); NaN(1,ad.nofig); pts(selected+1:size(pts,1),:)];
        pts(selected+1,cur_fig) = pts_new;
        shell = [shell(1:selected,:); NaN(1,ad.nofig); shell(selected+1:size(shell,1),:)];
        shell(selected+1,cur_fig) = shell_new;        
        selected = selected+1;
    else
        %new point
        points = [points; NaN(1,3*ad.nofig)];
        points(selected+1,3*(cur_fig-1)+1:3*cur_fig) = point;
        pts = [pts; NaN(1,ad.nofig)];
        pts(selected+1,cur_fig) = pts_new;
        shell = [shell; NaN(1,ad.nofig)];
        shell(selected+1,cur_fig) = shell_new;
        selected = selected+1;
    end
else
    %new point
    points = [points; NaN(1,3*ad.nofig)];
    points(selected+1,3*(cur_fig-1)+1:3*cur_fig) = point;
    pts = [pts; NaN(1,ad.nofig)];
    pts(selected+1,cur_fig) = pts_new;
    shell = [shell; NaN(1,ad.nofig)];
    shell(selected+1,cur_fig) = shell_new;    
    selected = selected+1;
end

%Store fields to ad
ad.pts(end+1,:) = {pts};
ad.shell(end+1,:) = {shell};
ad.points(end+1,:) = {points};
ad.selected(end+1) = selected;
ad.cur_fig(end+1) = cur_fig;

% Utilty function, delete a point
function ad = DeletePoint(ad)
%Load old data
selected = ad.selected(end);
if selected == 0
    return
end
points = ad.points{end,1};
pts = ad.pts{end,1};
shell = ad.shell{end,1};
cur_fig = ad.cur_fig(end);

disp('Before Delete');
disp('points');
disp(points);
disp('pts');
disp(pts);
disp('selected');
disp(selected);

points(selected,:) = NaN(1,3*ad.nofig);
pts(selected,:) = NaN*ones(1,ad.nofig);
shell(selected,:) = NaN*ones(1,ad.nofig);

if all(isnan(pts(selected,:)))
    points = points(setdiff(1:size(pts,1),selected),:);
    pts = pts(setdiff(1:size(pts,1),selected),:);
    shell = shell(setdiff(1:size(shell,1),selected),:);
    selected = selected-1;
    if ~selected
        if ~size(pts,1)
            disp('After Delete');
            disp('points');
            disp(points);
            disp('pts');
            disp(pts);
            disp('selected');
            disp(selected);
            ad.points(end+1,:) = {points};
            ad.pts(end+1,:) = {pts};
            ad.shell(end+1,:) = {shell};
            ad.cur_fig(end+1) = cur_fig;
            ad.selected(end+1) = selected;
            return
        else
            selected = size(pts,1);
        end
    end
end

%Look for new selection
if  ~all(isnan(pts(:,cur_fig))) % stay in current figure
    idx_tmp = ~isnan(pts(:,cur_fig));
    selected = find(idx_tmp(1:selected),1,'last');
    if ~size(selected,2)
        selected = find(idx_tmp,1,'last');
    end
elseif ~all(isnan(pts(selected,setdiff(1:ad.nofig,cur_fig))))
    idx_tmp = ~isnan(pts(selected,:));
    cur_fig = find(idx_tmp(1:cur_fig),1,'last');
    if ~size(cur_fig,2)
        cur_fig = find(idx_tmp,1,'last');
    end
elseif ~all(isnan(pts))
    [selected, cur_fig] = find(~isnan(pts),1,'first');
else
    selected = 0;
end

disp('After Delete');
disp('points');
disp(points);
disp('pts');
disp(pts);
disp('selected');
disp(selected);

%Store fields to ad
ad.points(end+1,:) = {points};
ad.pts(end+1,:) = {pts};
ad.shell(end+1,:) = {shell};
ad.cur_fig(end+1) = cur_fig;
ad.selected(end+1) = selected;


% Check how many points are placed
% function ad = CheckNumber(ad)



% Utilty function, moves a point
function ad  = MovePoint(ad, selection)
ad = DeletePoint(ad);
ad = AddPoint(ad,selection,false,false);

%Save the movement as one action
ad.points = ad.points(setdiff(1:end,end-1),:);
ad.pts = ad.pts(setdiff(1:end,end-1),:);
ad.shell = ad.shell(setdiff(1:end,end-1),:);
ad.selected = ad.selected(setdiff(1:end,end-1));

%Utility function, shows the current selection on the figure
function ad = RefreshSelection(ad)
%reset old points
selected = ad.selected(end);
cur_fig = ad.cur_fig(end);
if size(ad.pts,1)>=2
    pts = ad.pts{end-1,:};
    pts_tmp = pts(:);
    idx = ~isnan(pts_tmp);
    set(pts_tmp(idx),'Color',[0 0 0],'Visible','off','MarkerSize', 15);
end
%update new points
if size(ad.pts,1)~= 0
    pts = ad.pts{end,:};
    pts_tmp = pts(:);
    idx = ~isnan(pts_tmp);
    set(pts_tmp(idx),'Color',[0 0 0],'Visible','on','MarkerSize', 15);
    if selected~=0
        if ~isnan(pts(selected,cur_fig))
        	set(pts(selected,cur_fig),'Color',[0 1 0],'MarkerSize', 30);
        end
        pts_tmp = pts(selected,setdiff(1:ad.nofig,cur_fig));
        idx = ~isnan(pts_tmp);
        set(pts_tmp(idx),'Color',[0 0 1],'MarkerSize', 30);
    end
else
    return
end

%Utility function, toggles the visibility of the registered mesh
function ad = ToggleOverlay(ad)
    for i = setdiff(1:ad.nofig,2)
        obj_tmp = ad.obj_inv{i};
        if strcmp(get(obj_tmp,'Visible'),'on')
            set(obj_tmp,'Visible','off');
        else
            set(obj_tmp,'Visible','on');
        end
    end
   
%Utility funcitono, toggles the Facealphavalues
function ad = ToggleFaceAlpha(ad)
    for i = setdiff(1:ad.nofig,2)
        obj_tmp = ad.obj_inv{i};
        if get(ad.obj{2},'FaceAlpha') ~= 1
            set(obj_tmp,'FaceAlpha',0.4);
            set(ad.obj{2},'FaceAlpha',1);
        else
            set(obj_tmp,'FaceAlpha',1);
            set(ad.obj{2},'FaceAlpha',0.4);
        end
    end
    
%Utility function, undo the last action
function ad = UndoAction(ad)
    %Load and replace data
    ad.points = ad.points(1:end-1,:);
    
    pts = ad.pts{end,:};
    pts_tmp = pts(:);
    idx = ~isnan(pts_tmp);
    set(pts,'Visible','off');
    ad.pts = ad.pts(1:end-1,:);
    
    ad.selected = ad.selected(1:end-1);
    if size(ad.pts,1) ~= 0
        disp('After Undo');
        disp('conn');
        disp(ad.conn{end,1});
        disp('ln');
        disp(ad.ln{end,1});
        disp('directions');
        disp(ad.directions{end,1});
        disp('points');
        disp(ad.points{end,1});
        disp('pts');
        disp(ad.pts{end,1});
        disp('selected');
        disp(ad.selected(end));
    end
    ad = RefreshSelection(ad);

%Utility function, changes the current figure
function ad = ChangeFigure(ad)
    nofig = ad.nofig;
    fig = ad.figs;
    cur_fig = ad.cur_fig(end);
    %calculate next picture
    if nofig ==1
        return
    elseif cur_fig == nofig
        cur_fig = 1;
    else
        cur_fig = cur_fig+1;
    end
    figure(fig(cur_fig));
    set(fig(cur_fig), 'WindowButtonDownFcn', {@ButtonPress fig}, ...
         'KeyPressFcn', {@GUI_VisualisationKeyPress ad.ax(cur_fig) @KeyPress fig});
     if ~size(ad.pts,1)
        ad.cur_fig(end) = cur_fig;
     else
        ad.selected(end+1) = ad.selected(end);
        ad.points(end+1,:) = ad.points(end,:);
        ad.pts(end+1,:) = ad.pts(end,:);
        ad.shell(end+1,:) = ad.shell(end,:);
        ad.cur_fig(end+1) = cur_fig;
     end
     
