function normal = TRI_MeanNormal(varargin)
%TRI_MeanNormal  Calculate mean normal, weighted by triangle area.
%
%   Syntax:
%    normal = TRI_MeanNormal(F, V, normalise)
%    normal = TRI_MeanNormal(normals, normalise)
%
%   Input:
%    F:         M-by-3 array defining a surface triangle mesh. The rows
%               correspond to different triangles and the columns
%               correspond to the three vertices that make up each
%               triangle. The elements are row indices into V.
%    V:         N-by-3 array defining vertices. The rows correspond to
%               different vertices and the columns correspond to X-, Y- and
%               Z-coordinates. The elements are coordinate values.
%    normalise: Logical indicating whether or not the normals should be
%               normalised to unit length. Optional, defaults to true.
%    normals:   M-by-3 array defining non-normalized triangle normals. The
%               rows correspond to different triangles and the columns
%               correspond to X-, Y- and Z-components. The elements are
%               vector components. Optional, defaults to
%               TRI_Normals(F, V, false).
%
%   Output:
%    normal: 3-element row vector containing the mean triangle normal of
%            the mesh, weighted by triangle areas. The elements are X-, Y-
%            and Z-components.
%
%   Effect: This function will calculate the mean normal of the triangles
%   of a mesh, weighted by the areas of those triangles.
%
%   Dependencies: TRI_Normals.m
%
%   Known parents: Muscle_AlignRegion.m
%                  Muscle_CutByRegionBorder.m
%                  Muscle_SelectByRegionBorder.m
%                  Muscle_RegionCentroids.m


%Created on 22/01/2007 by Ward Bartels.
%WB, 16/05/2007: Speed improvements; added normalise argument.
%WB, 02/03/2009: Added input of normals.
%Stabile, fully functional.


%Check if normalise argument is present
lastarg = varargin{nargin};
normpresent = nargin>1 && islogical(lastarg);

%Calculate triangle normals if necessary <<TRI_Normals.m>>
if nargin-normpresent<2
    normals = varargin{1};
else
    normals = TRI_Normals(varargin{1}, varargin{2}, false);
end

%Calculate sum of normals, weighted by area
normal = sum(normals, 1);

%Normalise if requested
if ~normpresent || lastarg
    normal = normal/norm(normal, 2);
end