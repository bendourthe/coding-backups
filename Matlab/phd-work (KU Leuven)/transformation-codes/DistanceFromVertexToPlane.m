function [distances, normals] = DistanceFromVertexToPlane(V, beta, normals)
%DistanceFromVertexToPlane  Calculate distance from vertices to planes.
%wb20060123
%
%   Syntax:
%    [distances, normals] = DistanceFromVertexToPlane(V, beta)
%
%   Input:
%    V:    N-by-3 array containing vertex coordinates. Each row represents
%          a vertex; the first, second and third columns represent X-, Y-
%          and Z-coordinates respectively.
%    beta: 3*N-by-3 array which defines the planes by specifying vertex
%          coordinates. Each row represents a vertex that lies in a plane;
%          the first, second and third colums represent X-, Y- and Z-
%          coordinates respectively. The vertices in beta are grouped by
%          three; each group of vertices defines one plane.
%
%   Output:
%    distances: M-by-N array containing distances from each vertex to the
%               specified planes. Each row corresponds to the same row in
%               V; each column corresponds to a plane and to the same
%               column in normals.
%    normals:   3-by-N array containing the normals of the planes. The
%               first, second and third rows represent X-, Y- and Z-
%               coordinates respectively. Each column corresponds to a
%               plane and to the same column in distances.
%
%   Effect: This function will calculate the perpendicular distances from
%   the vertices in V to the planes defined by beta. It is possible to run
%   the calculation for multiple points in one go, by putting each point on
%   a different row in V. Also, multiple planes can be used by
%   concatenating their beta-arrays vertically.
%
%   Dependencies: TRI_Normals.m
%
%   Known parents: TRI_IntersectWithPlane.m
%                  TRI_IntersectWithBoundedPlane.m
%                  IntersectThreePlanes.m
%                  IntersectLineAndPlane.m
%                  GUI_IndicateMultiPlane.m

%Created on 19/12/2005 by Ward Bartels.
%WB, 23/01/2006: Distance to multiple planes can now be calculated in one
%                go.
%Stabile, fully functional.


%Calculate normals of planes <<TRI_Normals.m>>
if nargin <3
    normals = TRI_Normals(reshape(1:size(beta, 1), 3, []).', beta).';
end

%Calculate perpendicular distances
distances = V*normals-ones(size(V, 1), 1)*dot(beta(1:3:end,:).', normals, 1);