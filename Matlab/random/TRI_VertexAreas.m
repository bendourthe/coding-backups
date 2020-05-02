function areas = TRI_VertexAreas(F, arg2, double, numvert)
%TRI_VertexAreas  Calculate areas of triangles that surround each vertex.
%
%   Syntax:
%    areas = TRI_VertexAreas(F, V, double, numvert)
%    areas = TRI_VertexAreas(F, areas, double, numvert)
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
%    numvert: Number of vertices, equal to N. Optional, defaults to
%             size(V, 1) in the first syntax and max(F(:)) in the second.
%    areas:   M-element column vector defining double triangle areas. The
%             elements correspond to different triangles.
%
%   Output:
%    areas: N-element column vector containing the total area of the
%           triangles surrounding each vertex. The elements correspond to
%           different vertices.
%
%   Effect: For each vertex of the provided mesh, this function will
%   calculate the surrounding area by summing the areas of the neighbouring
%   triangles.
%
%   Dependencies: TRI_Areas.m
%
%   Known parents: Muscle_AlignRegion.m
%                  Femur_ShaftAxis.m
%                  Muscle_SelectDatabaseBone.m
%                  MC_SplitRegions.m
%                  Hip_JointCenter.m

%Created on 16/05/2007 by Ward Bartels.
%WB, 04/03/2009: Added input of areas and vertex count.
%Stabile, fully functional.


%Get double triangle areas and count vertices if necessary <<TRI_Areas.m>>
if size(arg2, 2)==1
    areas = arg2;
    if nargin<4
        numvert = max(max(F));
    end
else
    areas = TRI_Areas(F, arg2, true);
    if nargin<4
        numvert = size(arg2, 1);
    end
end

%Sum areas for each vertex
areas = accumarray(F(:), [areas; areas; areas], [numvert 1]);

%Divide by 2 unless specified otherwise
if nargin<3 || ~double
    areas = areas/2;
end