function vectors = NormaliseVectors(vectors, dim, norms)
%NormaliseVectors  Normalise vectors by their 2-norms.
%
%   Syntax:
%    vectors = NormaliseVectors(vectors, dim, norms)
%
%   Input:
%    vectors: N-dimensional array defining vectors. Dimension dim
%             corresponds to different components and all other dimensions
%             correspond to different vectors. The elements are vector
%             components.
%    dim:     Dimension on which the different vector components are
%             defined. Optional, defaults to ndims(vectors). Note that
%             ndims(vectors) is never less than 2.
%    norms:   N- or (N-1)-dimensional array containing the vectors' norms. 
%             Dimension dim is singleton and all other dimensions
%             correspond to different vectors. Optional, will be calculated
%             automatically if omitted.
%
%   Output:
%    vectors: N-dimensional array containing normal vectors. Dimension dim
%             corresponds to different components and all other dimensions
%             correspond to different vectors. The elements are vector
%             components.
%
%   Effect: This function will normalise the provided set of vectors by the
%   vectors' respective 2-norms, or alternatively by the provided set of
%   norms.
%
%   Dependencies: VectorNorms.m
%
%   Known parents: TRI_Normals.m
%                  TRI_Areas.m
%                  TRI_VertexNormals.m
%                  Muscle_CutByRegionBorder.m
%                  Muscle_SelectByRegionBorder.m
%                  Femur_Landmarks.m

%Created on 02/03/2009 by Ward Bartels.
%WB, 12/03/2009: Added dimension argument.
%Stabile, fully functional.


%Store size of vectors
vdims = ndims(vectors);

%Sum last dimension by default
if nargin<2
    dim = vdims;
end

%Calculate norms if not provided <<VectorNorms.m>>
if nargin<3
    norms = VectorNorms(vectors, dim);
end

%Normalise vectors
repdims = ones(1, vdims);
repdims(dim) = size(vectors, dim);
vectors = vectors./repmat(norms, repdims);