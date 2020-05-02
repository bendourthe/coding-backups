function varargout = GUI_ManipulateObjects(varargin)
%GUI_ManipulateObjects  Interactively manipulate plot objects.
%wb20070914
%
%   Syntax:
%    [tf, manfun] = GUI_ManipulateObjects(obj, rotcen, rotcon, trncon)
%    [tf, manfun] = GUI_ManipulateObjects(obj, manfun)
%    [R, t, T] = GUI_ManipulateObjects(h, action)
%
%   Input:
%    obj:    Column vector of handles to graphics objects that are to be
%            manipulated. Each object must be able to take an hgtransform
%            object as parent. All objects must be children of the same
%            axes. If obj is a single hgtransform object, it will be
%            manipulated directly instead of being assigned to a new
%            hgtransform object.
%    rotcen: 3-element row vector containing the X-, Y- and Z-coordinates
%            of the center of rotation. Optional, defaults to [0 0 0].
%    rotcon: N-by-3 array containing 0 to 3 rotation constraints. Each row
%            is one constraint axis, the three columns containn X-, Y- and
%            Z-coordinates. Rotation will be constrained so there is no
%            rotation around the axes specified in rotcon. Optional,
%            defaults to zeros(0, 3).
%    trncon: N-by-3 array containing 0 to 3 translation constraints. Each
%            row is one constraint axis, the three columns containn X-, Y-
%            and Z-coordinates. Translation will be constrained so there is
%            no translation along the axes specified in trncon. Optional,
%            defaults to zeros(0, 3).
%    manfun: N-by-1 cell array containing 1 to 3 function handles. The
%            three elements correspond to the left, right and middle (or
%            left + right) mouse button. Each handle references a function
%            that generates an incremental 4-by-4 transformation matrix
%            suitable to be multiplied with hgtransform's 'Matrix'
%            property, based on these arguments:
%             xys: 2-element row vector containing the position of the
%                  mouse pointer in pixels, when the mouse button was
%                  pressed down.
%             xy1: 2-element row vector containing the position of the
%                  mouse pointer in pixels, before motion occurred.
%             xy2: 2-element row vector containing the position of the
%                  mouse pointer in pixels, after motion occurred.
%             tf:  Handle to the hgtransform object being interactively
%                  manipulated.
%    h:      Handle to any object on the figure for which interactive
%            manipulation is currently active.
%    action: String defining an action. Can be set to:
%             'sample':  Returns the current transformation parameters
%                        without interrupting the interactive manipulation.
%             'suspend': Removes the callbacks set by the first or second
%                        syntax, but leaves the hgtransform object and the
%                        application data intact, in effect pausing the
%                        interactive manipulation.
%             'restore': Resumes after 'suspend'.
%             'end':     Restores the graphics objects to their situation
%                        before this function was run (without applying any
%                        transformation permanently).
%
%   Output:
%    tf:     Handle to the hgtransform object being interactively manipulated.
%    manfun: N-by-1 cell array containing 1 to 3 function handles. Same as
%            manfun in the input, except that it is automatically generated
%            in the first syntax.
%    R:      3-by-3 orthonormal rotation matrix.  Each row corresponds to
%            an input coordinate, each column corresponds to an output
%            coordinate.
%    t:      3-element translation row vector.
%    T:      4-by-4 homogeneous transformation matrix. Each row corresponds
%            to an input coordinate, each column corresponds to an output
%            coordinate - thus, it is the transpose of the hgtransform's
%            'Matrix' property.
%
%   Effect: This function will enable the user to modify the position and
%   scaling of a set of plot objects interactively. It achieves this by
%   assigning the provided objects to an hgtransform object; the 'Matrix'
%   property of this object can be set interactively.
%       Using the first syntax, the left button orbits the object if a
%   point sufficiently close to its center was clicked, and it rolls the
%   object if a point further away was clicked. The middle button
%   translates the object in the camera plane and the right button
%   translates perpendicular to the camera plane when the mouse is moved
%   vertically.
%       Using the second syntax, the three mouse buttons will activate up
%   to three functions referenced in manfun. This allows for more
%   flexibility in transforming the object (e.g. apply scaling, change
%   center of rotation), at the expense of being less straightforward.
%       The first and second syntax may be used while interactive
%   manipulation is already active. The mouse button behaviour will be
%   modified to comply with the new inputs. In order to achieve this, pass
%   the tf output from the first call to the second call's obj argument.
%       The third syntax may be used to retrieve the transformation matrix,
%   along with the principal rotation and translation (calculated by
%   DecomposeTransformation.m).
%
%   Dependencies: GUI_CoordinateExtents.m
%                 DecomposeTransformation.m
%                 GUI_ManipulateObjects.m (recursive)
%                 GUI_Pixel2Axes.m
%                 DistanceFromVertexToLine.m
%                 AxisRotation.m
%                 CenterTransformation.m
%
%   Known parents: GUI_ManipulateObjects.m (recursive)

%Created on 14/09/2007 by Ward Bartels.
%Stabile, fully functional.



%---------------%
% Main function %
%---------------%


%Enforce correct number of input arguments
error(nargchk(1, 4, nargin, 'struct'));

%Define pointer properties and callback properties
pointerprop = {'Pointer' 'PointerShapeCData' 'PointerShapeHotSpot'};
callbackprop = {'WindowButtonDownFcn' 'WindowButtonMotionFcn' 'WindowButtonUpFcn'};

%Check which syntax was used
if nargin<2 || ~ischar(varargin{2}) %first or second syntax - initialise
    
    %Get object, axes and figure handles
    obj = varargin{1};
    ax = ancestor(obj(1), 'axes');
    fig = ancestor(ax, 'figure');
    
    %Create hgtransform object, and put obj under it
    if numel(obj)==1 && strcmpi(get(obj, 'type'), 'hgtransform')
        tf = obj;
    else
        tf = hgtransform('Parent', ax);
        set(obj, 'Parent', tf);
    end
    
    %Check if manipulator functions are provided
    if nargin==2 && iscell(varargin{2}) %second syntax
        manfun = varargin{2};
    else %first syntax
        
        %Set defaults
        if nargin<2, rotcen = [0 0 0]; else rotcen = varargin{2}; end
        if nargin<3, rotcon = zeros(0, 3); else rotcon = varargin{3}; end
        if nargin<4, trncon = zeros(0, 3); else trncon = varargin{4}; end
        
        %Get size measure for objects <<GUI_CoordinateExtents.m>>
        extents = GUI_CoordinateExtents(ax);
        extents(isnan(extents)) = 0;
        center = mean(extents, 1);
        extent = mean(diff(extents, 1, 1))/2;
        if extent==0
            extent = get(ax, {'XLim' 'YLim' 'ZLim'});
            extent = mean(diff(vertcat(extent{:}), 1, 2));
        end
        radius = norm(center-rotcen)+extent;
        
        %Reverse transform rotcen and center (will be undone later)
        T = get(tf, 'Matrix');
        R = T(1:3,1:3).'; t = T(1:3,4).';
        rotcen = (rotcen-t)/R;
        center = (center-t)/R;
        
        %Calculate null space transformations for constraints
        n1 = null(rotcon);
        rotcon = n1*n1.';
        n2 = null(trncon);
        trncon = n2*n2.';
        
        %Create manipulator functions <Rotator> <Translator>
        manfun = {@(xys, xy1, xy2, tf) Rotator(xys, xy1, xy2, tf, ax, rotcen, radius, center, extent, rotcon);
                  @(ign1, xy1, xy2, ign2) Translator([0 0 xy1(2)-xy2(2)], ax, trncon);
                  @(ign1, xy1, xy2, ign2) Translator(xy2-xy1, ax, trncon)};
    end
    
    %Set pointer
    pointer = get(fig, pointerprop);
    setptr(fig, 'hand');
    
    %Store application data
    ad.tf = tf;
    ad.manfun = manfun;
    ad.pointer = pointer;
    ad.callbacks = get(fig, callbackprop);
    setappdata(fig, mfilename, ad);
    
    %Set callback <ButtonDown>
    set(fig, 'WindowButtonDownFcn', @ButtonDown);
    
    %Generate output
    varargout = {tf manfun};
    
else %third syntax - finalise
    
    %Get figure handle
    fig = ancestor(varargin{1}, 'figure');
    
    %Make sure interactive manipulation is active
    if ~isappdata(fig, mfilename)
        error([mfilename ':NotActive'], 'Interactive manipulation is not active.');
    end
    
    %Get application data
    ad = getappdata(fig, mfilename);
    
    %Read and decompose transformation matrix <<DecomposeTransformation.m>>
    T = get(ad.tf, 'Matrix').';
    [R, t] = DecomposeTransformation(T);
    
    %Check which action was requested
    switch varargin{2}
        case 'sample'
            
            %Do nothing
            
        case 'suspend'
            
            %Restore callbacks and reset pointer
            set(fig, callbackprop, ad.callbacks);
            set(fig, pointerprop, ad.pointer);
            
        case 'restore'
            
            %Re-initialise manipulation <<GUI_ManipulateObjects.m>>
            GUI_ManipulateObjects(ad.tf, ad.manfun);
            
        case 'end'
            
            %Restore callbacks and reset pointer
            set(fig, callbackprop, ad.callbacks);
            set(fig, pointerprop, ad.pointer);
            
            %Remove hgtransform object
            set(get(ad.tf, 'Children'), 'Parent', get(ad.tf, 'Parent'));
            delete(ad.tf);
            
            %Remove application data
            rmappdata(fig, mfilename);
            
        otherwise
            error([mfilename ':UnknownAction'], 'Unrecognised action.');
    end
    
    %Generate output
    varargout = {R t T};
end



%-----------%
% Callbacks %
%-----------%


%Callback, executed when a mouse button is pressed down
function ButtonDown(fig, evt)

%Simulate button-up event <ButtonUp>
ButtonUp(fig, evt);

%Set figure units to pixels and get current point
set(fig, 'Units', 'pixels');
xy1 = get(fig, 'CurrentPoint');

%Load application data
ad = getappdata(fig, mfilename);

%Store selection type
button = find(strcmpi(get(fig, 'SelectionType'), {'normal' 'alt' 'extend'}));
if numel(button)~=1 || button>numel(ad.manfun), return; end

%Set pointer
setptr(fig, 'closedhand');

%Update application data
ad.xys = xy1;
ad.xy1 = xy1;
ad.button = button;
ad.callbacks([2 3]) = get(fig, {'WindowButtonMotionFcn' 'WindowButtonUpFcn'});
setappdata(fig, mfilename, ad);

%Set button motion and release callback <ButtonMotion> <ButtonUp>
set(fig, 'WindowButtonMotionFcn', @ButtonMotion);
set(fig, 'WindowButtonUpFcn', @ButtonUp);



%Callback, executed when the mouse cursor is moved within the figure
function ButtonMotion(fig, evt)

%Set figure units to pixels and get current point
set(fig, 'Units', 'pixels');
xy2 = get(fig, 'CurrentPoint');

%Load application data
ad = getappdata(fig, mfilename);

%Run manipulator function and update transformation
T = ad.manfun{ad.button}(ad.xys, ad.xy1, xy2, ad.tf);
set(ad.tf, 'Matrix', T*get(ad.tf, 'Matrix'));

%Update application data
ad.xy1 = xy2;
setappdata(fig, mfilename, ad);



%Callback, executed when the mouse button is released
function ButtonUp(fig, evt)

%Set pointer
setptr(fig, 'hand');

%Get application data
ad = getappdata(fig, mfilename);

%Restore figure button motion and button release callback
set(fig, {'WindowButtonMotionFcn' 'WindowButtonUpFcn'}, ad.callbacks([2 3]));



%-------------------%
% Utility functions %
%-------------------%


%Utility function, creates rotation matrix from mouse motion
function T = Rotator(xys, xy1, xy2, tf, ax, rotcen, radius, center, extent, constr)

%Update position of rotcen (may be influenced by translation)
T_old = get(tf, 'Matrix').';
R = T_old(1:3,1:3);
t = T_old(4,1:3);
rotcen = rotcen*R+t;
center = center*R+t;

%Calculate selected points in axes coordinates <<GUI_Pixel2Axes.m>>
points = GUI_Pixel2Axes(ax, [xys; xy1; xy2], false);

%Calculate displacement and viewing plane normal
displ = points(3,:)-points(2,:);
viewdir = campos(ax)-camtarget(ax);
viewdir = viewdir/norm(viewdir);

%Check if clicked point is outside range <<DistanceFromVertexToLine.m>>
if DistanceFromVertexToLine(points(1,:), center, viewdir)>extent
    
    %Determine rotation axis and angle for "rolling" motion
    rotax = viewdir;
    moment = cross(viewdir, points(2,:)-rotcen);
    rotang = atan((moment*displ.')/norm(moment)^2);
    
else
    
    %Determine rotation axis and angle for "orbiting" motion
    rotax = cross(viewdir, displ);
    rotang = atan(norm(displ)/radius);
end

%Fail-safe: return unity transformation if rotation is too small
n_rotax = norm(rotax);
if n_rotax<eps
    T = eye(4);
    return
end

%Remove displacement components parallel to constr
rotax = rotax/n_rotax*rotang;
rotax = rotax*constr;
rotang = norm(rotax);

%Fail-safe: return unity transformation if rotang is too small
if rotang<eps
    T = eye(4);
    return
end

%Create transformation matrix for hgtransform <<AxisRotation.m>>
%                                             <<CenterTransformation.m>>
R = AxisRotation(rotax, rotang);
t = CenterTransformation(R, rotcen);
T = [R.' t.'; 0 0 0 1];



%Utility function, creates translation matrix from mouse motion
function T = Translator(displ, ax, constr)

%Transform displacement into axes coordinates <<GUI_Pixel2Axes.m>>
displ = GUI_Pixel2Axes(ax, displ, true);

%Remove displacement components parallel to constraints
displ = displ*constr;

%Create transformation matrix for hgtransform
T = [eye(3) displ.'; 0 0 0 1];