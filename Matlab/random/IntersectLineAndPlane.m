function vertices = IntersectLineAndPlane(alpha, beta)
%IntersectLineAndPlane  Calculate intersection points of lines and planes.
%wb20060619
%
%   Syntax:
%    vertices = IntersectLineAndPlane(alpha, beta)
%
%   Input:
%    alpha: 2*N-by-3 array which defines the lines by specifying vertex
%           coordinates. Each row represents a vertex that lies on a line;
%           the first, second and third colums represent X-, Y- and Z-
%           coordinates respectively. The vertices in alpha are grouped by
%           two; each group of vertices defines one line.
%    beta:  3*N-by-3 array which defines the planes by specifying vertex
%           coordinates. Each row represents a vertex that lies in a plane;
%           the first, second and third colums represent X-, Y- and Z-
%           coordinates respectively. The vertices in beta are grouped by
%           three; each group of vertices defines one plane.
%
%   Ouptut:
%    vertices: M-by-3-by-N array containing vertex coordinates. The first
%              dimension corresponds to different lines in alpha, the
%              second dimension corresponds to X-, Y- and Z-coordinates and
%              the third dimension corresponds to different planes in beta.
%
%   Effect: This function will calculate the intersection points of a set
%   of planes and a set of lines. Multiple lines and planes can be defined
%   by concatenating their alpha- or beta-arrays vertically.
%
%   Dependencies: DistanceFromVertexToPlane.m
%
%   Known parents: GUI_IndicateMultiPlane.m
%                  TRI_IntersectWithLine.m

%Created on 13/03/2006 by Ward Bartels.
%WB, 19/06/2006: Multiple planes can now be intersected.
%Stabile, fully functional.


%Get points from alpha
a1 = alpha(1:2:end-1,:);
a2 = alpha(2:2:end,:);

%Calculate distances from line points to plane
d1 = DistanceFromVertexToPlane(a1, beta);
d2 = DistanceFromVertexToPlane(a2, beta);

%Calculate intersection points
%Equivalent (single plane): vertices = a1+(a2-a1).*((d1./(d1-d2))*[1 1 1])
rep = [1 1 size(d1, 2)];
vertices = repmat(a1, rep)+repmat(a2-a1, rep).*repmat(permute(d1./(d1-d2), [1 3 2]), [1 3 1]);