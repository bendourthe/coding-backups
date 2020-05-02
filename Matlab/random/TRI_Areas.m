function [areas, normals] = TRI_Areas(F, V, double)
%TRI_Areas  Calculate triangle areas.
%
%   Syntax:
%    [areas, normals] = TRI_Areas(F, V, double)
%
%   Input:
%    F:       M-by-3 array defining a surface triangle mesh. The rows
%             correspond to different triangles and the columns correspond
%             to the three vertices that make up each triangle. The
%             elements are row indices into V.
%    V:       N-by-3 array defining vertices. The rows correspond to
%             different vertices and the columns correspond to X-, Y- and
%             Z-coordinates. The elements are coordinate values.
%    double:  Logical indicating whether or not the double of the areas
%             should be returned. This is slightly faster, and may be
%             useful when the areas act as weights. Optional, defaults to
%             false.
%
%   Output:
%    areas:   M-element column vector containing triangle areas. The
%             elements correspond to different triangles.
%    normals: M-by-3 array containing normalized triangle normals. The rows
%             correspond to different triangles and the columns correspond
%             to X-, Y- and Z-components. The elements are vector
%             components.
%
%   Effect: This function will calculate the area of each triangle in the
%   provided mesh. The normalized triangle normals are returned if
%   requested.
%
%   Dependencies: TRI_Normals.m
%                 VectorNorms.m
%                 NormaliseVectors.m
%
%   Known parents: TRI_VertexAreas.m
%                  KIN_SurfaceRotationCenter.m
%                  Femur_Landmarks.m
%                  KIN_JointSurfaceDistance.m
%                  Muscle_SplitRegion.m

%Created on 16/05/2007 by Ward Bartels.
%WB, 18/12/2007: Added output of normals.
%Stabile, fully functional.


%Calculate normals (not normalized) <<TRI_Normals.m>>
normals = TRI_Normals(F, V, false);

%Calculate double areas from normals <<VectorNorms.m>>
areas = VectorNorms(normals);

%Normalise normals if requested <<NormaliseVectors.m>>
if nargout>=2
    normals = NormaliseVectors(normals, 2, areas);
end

%Divide by 2 unless specified otherwise
if nargin<3 || ~double
    areas = areas./2;
end