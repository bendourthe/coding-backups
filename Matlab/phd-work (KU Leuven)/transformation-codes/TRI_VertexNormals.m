function normals = TRI_VertexNormals(F, V, normalise, normals)
%TRI_VertexNormals  Calculate vertex normals, weighted by triangle area.
%
%   Syntax:
%    normals = TRI_VertexNormals(F, V, normalise, normals)
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
%    normals: N-by-3 array containing vertex normals. The rows correspond
%             to different vertices and the columns correspond to X-, Y-
%             and Z-components. The elements are vector components.
%
%   Effect: This function will calculate vertex normals for the provided
%   mesh. For polyhedra inscribed in a sphere, the resulting vertex normals
%   will be analytically correct. This function will only run properly if
%   all redundancy has been removed from the mesh.
%
%   References: Max, N. Weights for computing vertex normals from facet
%               normals. J Graph Tools 4(2), 1999, 1-6.
%
%   Dependencies: TRI_Normals.m
%                 NormaliseVectors.m
%
%   Known parents: Muscle_AlignRegion.m
%                  Muscle_CutByRegionBorder.m
%                  Femur_ShaftAxis.m
%                  SIMM_WriteModel.m
%                  Hip_JointCenter.m

%Created on 06/12/2006 by Ward Bartels.
%WB, 08/05/2007: Speed improvements.
%WB, 16/05/2007: Added normalise argument.
%WB, 03/03/2009: Changed algorithm to Max's published method.
%Stabile, fully functional.


%Calculate normals if necessary <<TRI_Normals.m>>
if nargin<4
    normals = TRI_Normals(F, V, false);
end

%Calculate squared edge lengths
vertices = reshape(V(F,:), size(F, 1), 3, 3);
selen =  sum((vertices(:,[3 1 2],:)-vertices(:,[2 3 1],:)).^2, 3);

%Calculate weights for the normals
weights = 1./(selen(:,[3 1 2]).*selen(:,[2 3 1]));

%Accumulate weighted triangle normal components for all vertices
F = F(:);
one_array = ones(size(F));
weights = weights(:);
normals = [normals; normals; normals].*[weights weights weights];
normals = accumarray([F one_array; F one_array*2; F one_array*3], normals(:), size(V));

%Normalise normals if requested <<NormaliseVectors.m>>
if nargin<3 || normalise
    normals = NormaliseVectors(normals);
end


% %Original code (wb20090302), begin:
% 
% %Set default for normalise
% if nargin<3, normalise = true; end
% 
% %Calculate normals (non-normalized, already weighted) <<TRI_Normals.m>>
% normals = TRI_Normals(F, V, false);
% 
% %Accumulate triangle normal components for all vertices
% F = F(:);
% one_array = ones(size(F));
% normals = [normals; normals; normals];
% normals = accumarray([F one_array; F one_array*2; F one_array*3], normals(:));
% 
% %Calculate projected areas if requested or needed <<VectorNorms.m>>
% if normalise || nargout>=2
%     areas = VectorNorms(normals);
% end
% 
% %Normalise if requested
% if normalise
%     normals = normals./(areas*[1 1 1]);
% end
% 
% %Divide parallellogram areas by 2
% if nargout>=2
%     areas = areas/2;
% end
% 
% %Original code (wb20090302), end.