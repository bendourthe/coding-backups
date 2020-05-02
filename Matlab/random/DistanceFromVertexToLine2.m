function val = DistanceFromVertexToLine2(V, orig, vect)
%DistanceFromVertexToLine  Calculate distance from vertices to lines.
%
%   Syntax:
%    distances = DistanceFromVertexToLine(V, orig, vect)
%
%   Input:
%    V:    1-by-3 array defining vertex. The rows correspond to different
%          vertices and the columns correspond to X-, Y- and Z-coordinates.
%          The elements are coordinate values.
%    orig: N-by-3 array defining vertices that are located on the lines.
%          The rows correspond to different lines and the columns
%          correspond to X-, Y- and Z-coordinates. The elements are
%          coordinate values.
%    vect: N-by-3 array defining normal vectors that indicate the lines'
%          directions. The rows correspond to different lines and the
%          columns correspond to X-, Y- and Z-components. The elements are
%          vector components.
%
%   Output:
%    val: sum of squares of distances from V to the lines
%   
%   Effect: This function will calculate the perpendicular distances from
%   the vertices in V to the lines defined by orig and vect. Distances are
%   calculated for each combination of a vertex and a line.
%
%   Dependencies: VectorNorms.m
%
%   Known parents: GUI_ManipulateObjects.m

%Created on 23/01/2006 by Ward Bartels.
%WB, 09/03/2009: Changed line input to origin and vector representation.
%Stabile, fully functional.


% %Rearrange and replicate orig and vect to conform with vertices
% orig = repmat(reshape(orig, [1 size(orig, 1) 3]), [size(V, 1) 1 1]);
% vect = repmat(reshape(vect, [1 size(vect, 1) 3]), [size(V, 1) 1 1]);
% 
% %Rearrange and replicate vertices to conform with orig and vect
% V = repmat(reshape(V, [size(V, 1) 1 3]), [1 size(orig, 2) 1]);

V = repmat(V, size(orig, 1), 1);

%Calculate distances <<VectorNorms.m>>
distances = realsqrt(sum(realpow((cross(V-orig, vect, 2)), 2), 2));
val = sum(realpow(distances, 2), 1);