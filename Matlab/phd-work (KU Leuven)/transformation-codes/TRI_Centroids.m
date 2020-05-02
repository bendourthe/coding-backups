function centroids = TRI_Centroids(F, V)
%TRI_Centroids  Determine centroids of mesh triangles or mesh tetrahedrons.
%
%   Syntax:
%    centroids = TRI_Centroids(F, V)
%
%   Input:
%    F: M-by-3/4 array defining a surface triangle mesh. The rows correspond
%       to different triangles and the columns correspond to the three
%       vertices that make up each triangle. The elements are row indices
%       into V.
%    V: N-by-3 array defining vertices. The rows correspond to different
%       vertices and the columns correspond to X-, Y- and Z-coordinates.
%       The elements are coordinate values.
%
%   Output:
%    centroids: M-by-3 array containing triangle centroids. The rows
%               correspond to different triangles and the columns
%               correspond to X-, Y- and Z-coordinates. The elements are
%               coordinate values.
%
%   Effect: This function calculates the centroids of all triangles. These
%   are calculated as the mean coordinates of the triangles' corner points.
%
%   Dependencies: none
%
%   Known parents: Muscle_SelectByRegionBorder.m
%                  TRI_DetermineTriangleConnections.m
%                  KIN_SurfaceRotationCenter.m
%                  TRI_FlipNormalsToConvex.m
%                  Femur_Landmarks.m
%                  TRI_CutAlongContour.m
%                  Muscle_SplitRegion.m

%Created on 25/04/2007 by Ward Bartels.
%Adjusted on 10/05/2011 by Steven Walscharts (expand to tetrahedrons)
%Stabile, fully functional.


%Calculate centroids
centroids = permute(mean(reshape(V(F,:), size(F, 1), size(F, 2), 3), 2), [1 3 2]);