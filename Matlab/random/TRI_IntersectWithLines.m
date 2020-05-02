function [vertices, triangles, param] = TRI_IntersectWithLines(F, V, orig, vect, sortfun, sparse_out)
%TRI_IntersectWithLines  Intersect mesh and lines.
%wb20070301
%
%   Syntax:
%    [vertices, triangles, param] = ...
%        TRI_IntersectWithLines(F, V, orig, vect, sortfun, sparse_out)
%
%   Input:
%    F:          N-by-3 array containing indices into V. Each row
%                represents a triangle, each element is a link to a vertex
%                in V.
%    V:          N-by-3 array containing vertex coordinates. Each row
%                represents a vertex; the first, second and third columns
%                represent X-, Y- and Z-coordinates respectively.
%    orig:       N-by-3 array containing ray origin coordinates. Each row
%                corresponds to a ray; the three columns contain X-, Y- and
%                Z- coordinates.
%    vect:       N-by-3 array containing ray direction vector components. Each
%                row corresponds to a ray; the three columns contain X-, Y-
%                and Z-components.
%    sortfun:    Function handle to a function that takes param, ray,
%                triangles and vertices as arguments, and returns a column
%                vector. The ray argument is a column vector containing row
%                indices into orig and vect. Optional, defaults to
%                @(p, r, t, v) -abs(p(:,1)).
%    sparse_out: Logical indicating whether or not the output arrays should
%                be sparse. Optional, defaults to false.
%
%   Output:
%    vertices:  N-by-3 array containing vertex coordinates. Each row
%               represents the intersection of one line with the mesh; the
%               first, second and third columns represent X-, Y- and Z-
%               coordinates respectively.
%    triangles: Column vector containing the intersected triangles. Each
%               element is a row index into F.
%    param:     N-by-3 array containing parametric coordinates of
%               intersection vertices. Each row corresponds to the
%               intersection of one line with one triangle. The first
%               column contains parametric ray coordinates, the second and
%               third column contain barycentric triangle coordinates.
%
%   Effect: This function will intersect triangles in F with the lines
%   defined by orig and vect. For lines that intersect the mesh in multiple
%   points, only the intersection point yielding the highest output through
%   sortfun will be kept. Intersections for which the corresponding element
%   in sortfun's output is NaN are deleted. When a line has no more
%   intersections, the corresponding rows in vertices and param will
%   contain NaN if sparse_out is false and 0 if it is true. The
%   corresponding elements in triangles will contain 0.
%
%   Dependencies: TRI_IntersectWithRays.m
%                 StackEqualElementIndices.m
%
%   Known parents: TRI_ShellSeparationDistance.m
%                  TRI_IntersectWithVectors.m
%                  Muscle_CutByRegionBorder.m

%Created on 13/07/2006 by Ward Bartels.
%WB, 01/03/2007: Rewrite based on TRI_IntersectWithRays.m.
%Stabile, fully functional.


%Set default inputs
if nargin<5, sortfun = @(p, r, t, v) -abs(p(:,1)); end
if nargin<6, sparse_out = false; end

%Intersect mesh with rays <<TRI_IntersectWithRays.m>>
[param, ray, triangles, vertices] = TRI_IntersectWithRays(F, V, orig, vect);

%Calculate ranking number
rank = sortfun(param, ray, triangles, vertices);

%Find highest rank per ray <<StackEqualElementIndices.m>>
[rank, ind] = sort(rank, 1, 'ascend');
ind(isnan(rank),:) = []; %Remove ranks equal to NaN
ind = ind(StackEqualElementIndices(ray(ind), 0, 1));

%Keep only highest ranking intersections per ray
vertices = vertices(ind,:);
triangles = triangles(ind,:);
ray = ray(ind,:);
param = param(ind,:);

%Expand ray index vector and create coordinate indices
ind_ray = ray(:,[1 1 1]);
ind_coord = repmat([1 2 3], size(ray, 1), 1);

%Check if sparse output was requested and create output arrays
if sparse_out
    vertices = sparse(ind_ray, ind_coord, vertices, size(orig, 1), 3);
    triangles = sparse(ray, ones(size(ray)), triangles, size(orig, 1), 1);
    param = sparse(ind_ray, ind_coord, param, size(orig, 1), 3);
else
    ind_ray = ind_ray(:);
    ind_coord = ind_coord(:);
    vertices = accumarray([ind_ray ind_coord], vertices(:), [size(orig, 1) 3], [], NaN);
    triangles = accumarray(ray, triangles, [size(orig, 1) 1]);
    param = accumarray([ind_ray ind_coord], param(:), [size(orig, 1) 3], [], NaN);
end


% %Original code (wb20070301), begin:
% 
% %Initialise output variables
% N = size(alpha, 1)/2;
% vertices = NaN(N, 3);
% if nargout>=2
%     triangles = zeros(N, 1);
% end
% 
% %Get coordinates of triangle corner points
% p1 = [V(F(:,1),:) ones(size(F, 1), 1)];
% p2 = [V(F(:,2),:) ones(size(F, 1), 1)];
% p3 = [V(F(:,3),:) ones(size(F, 1), 1)];
% 
% %Loop over all lines
% for ind = 1:N
%     
%     %Assemble new X-, Y- and Z-axes for projection
%     zax = alpha(2*ind,:)-alpha(2*ind-1,:);
%     [ignoble, jnd] = min(abs(zax));
%     xax(jnd) = 0;
%     xax(mod(jnd, 3)+1) = zax(mod(jnd+1, 3)+1);
%     xax(mod(jnd+1, 3)+1) = -zax(mod(jnd, 3)+1);
%     yax(jnd) = -zax(mod(jnd, 3)+1)^2-zax(mod(jnd+1, 3)+1)^2;
%     yax(mod(jnd, 3)+1) = zax(jnd)*zax(mod(jnd, 3)+1);
%     yax(mod(jnd+1, 3)+1) = zax(jnd)*zax(mod(jnd+1, 3)+1);
%     
%     %Project points
%     tform = inv([xax; yax; zax]);
%     tform(4,:) = -alpha(2*ind-1,:)*tform; %Origin at first line point
%     tf = tform(:,1:2);
%     t1 = p1*tf;
%     t2 = p2*tf;
%     t3 = p3*tf;
%     
%     %Get rotation directions of triangle side vectors around origin
%     rotdir = [t2(:,1).*t3(:,2)-t3(:,1).*t2(:,2) ...
%               t3(:,1).*t1(:,2)-t1(:,1).*t3(:,2) ...
%               t1(:,1).*t2(:,2)-t2(:,1).*t1(:,2)];
%     
%     %A triangle is intersected when rotdir contains negative or positive
%     %numbers, but not both, on one row
%     intersected = xor(any(rotdir>eps, 2), any(rotdir<-eps, 2));
%     
%     %Calculate intersection points' Z-coordinates
%     V1 = [t1(intersected,1) t2(intersected,1) t3(intersected,1)];              %X-coordinates in new system
%     V2 = [t1(intersected,2) t2(intersected,2) t3(intersected,2)];              %Y-coordinates in new system
%     V3 = [p1(intersected,:); p2(intersected,:); p3(intersected,:)]*tform(:,3); %Z-coordinates in new system
%     V3 = reshape(V3, [], 3);
%     z = (V1(:,1).*(V2(:,2).*V3(:,3)-V2(:,3).*V3(:,2))+...
%          V1(:,2).*(V2(:,3).*V3(:,1)-V2(:,1).*V3(:,3))+...
%          V1(:,3).*(V2(:,1).*V3(:,2)-V2(:,2).*V3(:,1)))./...
%         (V1(:,1).*(V2(:,2)-V2(:,3))+...
%          V1(:,2).*(V2(:,3)-V2(:,1))+...
%          V1(:,3).*(V2(:,1)-V2(:,2)));
%     
%     %Keep intersection point closest to first line point (origin)
%     [ignoble, jnd] = min(abs(z));
%     if ~isempty(jnd)
%         vertices(ind,:) = alpha(2*ind-1,:)+z(jnd)*zax;
%     end
%     
%     %Determine index of selected triangle, if necessary
%     if nargout>=2
%         triangles_intersected = find(intersected);
%         if isempty(triangles_intersected)
%             triangles(ind, 1) = 0;
%         else
%             triangles(ind, 1) = triangles_intersected(jnd);
%         end
%     end
% end
% 
% %Original code (wb20070301), end.