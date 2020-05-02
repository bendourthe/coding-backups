function [loc, trans] = KIN_ScrewAxis(rot, trans, loc)
%KIN_ScrewAxis  Locate screw axes from translation and rotation vectors.
%
%   Syntax:
%    loc = KIN_ScrewAxis(rot, trans, loc)
%
%   Input:
%    rot:   N-by-3 array defining rotation velocity vectors. The rows
%           correspond to different screw axes and the columns correspond
%           to X-, Y- and Z-components. The elements are vector components
%           expressed in radians per time unit.
%    trans: N-by-3 array defining translation velocity vectors. The rows
%           correspond to different screw axes and the columns correspond
%           to X-, Y- and Z-components. The elements are vector components
%           expressed in distance units per time unit.
%    loc:   N-by-3 array defining the points for which trans is the
%           translation velocity. The rows correspond to different screw
%           axes and the columns correspond to X-, Y- and Z-coordinates.
%           The elements are coordinate values. Optional, defaults to
%           zeros(size(rot, 1), 3).
%
%   Output:
%    loc:   N-by-3 array containing points on the screw axes. The rows
%           correspond to different screw axes and the columns correspond
%           to X-, Y- and Z-coordinates. The elements are coordinate
%           values.
%
%   Effect: This function will locate the instantaneous screw axes of the
%   motions of a set of rigid bodies, based on their rotation vectors and
%   the velocity of one point of each body.
%
%   Dependencies: none
%
%   Known parents: KIN_DirectedSpringRotation.m

%Created on 29/01/2008 by Ward Bartels.
%WB, 10/03/2009: Removed output of trans.
%Stabile, fully functional.


%Calculate vector towards point on screw axis
displ = cross(rot, trans, 2)./(sum(rot.^2, 2)*[1 1 1]);

%Displace loc
if nargin<3
    loc = displ;
else
    loc = loc+displ;
end