function distances = SquaredDistanceFromVertexToVertex(V1, V2, paired)
%SquaredDistanceFromVertexToVertex  Calculate squared vertex distances.
%
%   Syntax:
%    distances = SquaredDistanceFromVertexToVertex(V1, V2, paired)
%
%   Input:
%    V1:     N-by-3 array containing vertex coordinates. Each row
%            represents a vertex; the first, second and third columns
%            represent X-, Y- and Z-coordinates respectively.
%    V2:     N-by-3 array containing vertex coordinates. Each row
%            represents a vertex; the first, second and third columns
%            represent X-, Y- and Z-coordinates respectively. Optional,
%            defaults to V1.
%    paired: Logical indicating whether or not the points in V1 and V2 are
%            paired. If set to true, V1 and V2 must have the same sizes and
%            for each point in V1, only the squared distance to the point
%            in V2 on the same row will be calculated. Optional, defaults
%            to false.
%
%   Output:
%    distances: M-by-N array containing squared distances between points.
%               If paired is set to false, each row corresponds to a row in
%               V1 and each column corresponds to a row in V2; otherwise,
%               distances is a column vector and each element corresponds
%               to a row in V1 and V2.
%
%   Effect: This function will calculate the squared distances between
%   points in point clouds. If one point cloud is provided, it will
%   calculate all distances between points in that cloud. If two clouds are
%   provided and paired is set to false, it will calculate the distances
%   between every possible combination of points from both clouds. If
%   paired is set to true, the point clouds are considered to consist of
%   matched points, and only distances between matched points will be
%   calculated.
%
%   Dependencies: none
%
%   Known parents: ICP_RegisterMatchedDatasets.m
%                  ICP_MatchDatasets.m
%                  Femur_Landmarks.m

%Created on 28/11/2006 by Ward Bartels.
%Stabile, fully functional.


%Distinguish between possible inputs
if nargin<2 %One point cloud provided
    
    %Calculate distance for each possible pair of points from the cloud
    products = V1*V1.';
    distances = diag(products);
    distances = distances(:,ones(1, size(V1, 1)));
    distances = distances+distances.'-2*products;

elseif nargin>=3 && paired %Two paired point clouds provided
    
    %Calculate distance for each pair of points
    distances = sum((V2-V1).^2, 2);
    
elseif size(V1, 1)~=1 && size(V2, 1)~=1 %Two point clouds provided
    
    %Transpose V2
    V2 = V2.';
    
    %Calculate distance for each possible pair of points from both clouds
    squares1 = sum(V1.^2, 2);
    squares2 = sum(V2.^2, 1);
    distances = squares1(:,ones(size(squares2)))+...
                squares2(ones(size(squares1)),:)-...
                2*(V1*V2);
    
else %Point cloud and single point provided
    
    %Calculate distance from point to every point in the cloud
    if size(V2, 1)==1
        distances = sum(V1.^2, 2)+V2*V2.'-2*(V1*V2.');
    else
        distances = (sum(V2.^2, 2)+V1*V1.'-2*(V2*V1.')).';
    end
end