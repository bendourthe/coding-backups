function [minproj,mindist] = ClosestEdgePoint(point,vertices)
%ClosestEdgePoint Calculates the closest point on the edge of a contour
%
%   Syntax:
%    [minproj,mindist] = ClosestEdgePoint(point,vertices)
%
%   Input:
%    point:     3-column vector with the location of a point on the
%               triangle
%    vertices:  N-by-3 array containing vertex coordinates of the corners
%               of the contour
%
%   Output:
%    minproj:     3-column vector with the location of the point on the
%                 contour closest to the given point
%    mindist:     the distance of the point to the closest edge of the
%                 contour

sze = size(vertices,1);
dist = zeros(sze,1);
proj = zeros(sze,3);
orig = vertices;
vect = circshift(vertices,[-1 0]) - vertices;
vect = vect./(VectorNorms(vect)*ones(1,3)); %%%%%%%%%%%%%% VectorNorms(vect,2)
for i = 1:sze
    P = vect(i,:)'*vect(i,:);
    proj(i,:) = (orig(i,:)' + P*(point' - orig(i,:)'))';
    dist(i) = norm(point-proj(i,:));
end
[mindist, ind] = min(dist);
minproj = proj(ind,:);
end
