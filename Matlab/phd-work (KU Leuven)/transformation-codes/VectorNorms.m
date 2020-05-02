function norms = VectorNorms(vectors)
%VectorNorms  Calculate vector 2-norms (Euclidean length or distance).
%wb20070516
%
%   Syntax:
%    norms = VectorNorms(vectors)
%
%   Input:
%    vectors: N-by-M array containing vectors. Each row represents one
%             vector; the columns contain coordinate data in different
%             dimensions.
%
%   Output:
%    norms: Column vector containing the vectors' norms. Each element
%           corresponds to a row in vectors.
%
%   Effect: This function will calculate the 2-norms (regular Euclidean
%   length or distance) of the provided set of vectors.
%
%   Dependencies: none
%
%   Known parents: TRI_Normals.m
%                  TRI_VertexNormals.m
%                  TRI_Areas.m
%                  ICP_Stability.m
%                  GUI_Pixel2AxesTransform.m
%                  Muscle_AlignRegion.m

%Created on 16/05/2007 by Ward Bartels.
%Stabile, fully functional.


%Calculate the vector norms, using time-efficient functions.
norms = realsqrt(sum(realpow(vectors, 2), 2));