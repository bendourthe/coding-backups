function [V_new, contour, triangles, edges_cut, failed_test, distances] = TRI_IntersectWithBoundedPlane(F, V, beta, beta_bound)
%TRI_IntersectWithBoundedPlane  Intersect mesh with bounded plane.
%
%   Syntax:
%    [V_new, contour, triangles, edges_cut, failed_test, distances] = ...
%        TRI_IntersectWithBoundedPlane(F, V, beta, beta_bound)
%
%   Input:
%    F:          N-by-3 array containing indices into V. Each row
%                represents a triangle, each element is a link to a vertex
%                in V.
%    V:          N-by-3 array containing vertex coordinates. Each row
%                represents a vertex; the first, second and third columns
%                represent X-, Y- and Z-coordinates respectively.
%    beta:       3-by-3 array which defines the plane by specifying vertex
%                coordinates. Each row represents a vertex that lies in the
%                plane; the first, second and third colums represent X-, Y-
%                and Z-coordinates respectively.
%    beta_bound: 3*N-by-3 array which defines the boundary planes by
%                specifying vertex coordinates. Each row represents a
%                vertex that lies in a plane; the first, second and third
%                colums represent X-, Y- and Z- coordinates respectively.
%                The vertices in beta_bound are grouped by three; each
%                group of vertices defines one plane.
%
%   Output:
%    V_new:       N-by-3 array containing vertex coordinates of inter-
%                 section points. Each row represents a vertex; the first,
%                 second and third columns represent X-, Y- and Z-
%                 coordinates respectively.
%    contour:     N-by-2 array containing the intersection line. Each row
%                 represents an edge and corresponds to the same row in
%                 triangles; each element is a row index into V_new.
%    triangles:   Column vector containing the intersected triangles. Each
%                 element is a row index into F.
%    edges_cut:   N-by-2 array indicating which edge of a triangle is
%                 intersected by the plane. Each row corresponds to the
%                 same row in triangles and each column corresponds to the
%                 same column in contour. A 1 indicates the edge connecting
%                 columns 1 and 2 in F, a 2 indicates the edge connecting
%                 columns 2 and 3, and so on. For triangles intersected in
%                 only one edge, edges_cut will contain a zero in the
%                 second column; both columns will contain a zero when
%                 intersecting inside a triangle.
%    failed_test: N-by-2 array indicating which boundary test failed. Each
%                 row corresponds to the same row in triangles and each
%                 column corresponds to the same column in contour; the
%                 elements indicate which boundary plane truncated the
%                 corresponding contour edge. A zero indicates no
%                 truncation occurred.
%    distances:   Column vector containing distances from each vertex to
%                 the specified plane. Each row corresponds to the same row
%                 in V.
%
%   Effect: This function will extract all edges from the mesh and
%   determine which edges intersect the plane. For these edges, the
%   intersection points are calculated. The intersection line thus obtained
%   is truncated where it crosses one of the boundary planes defined by
%   beta_bound. The coordinates of the intersection points are returned
%   along with indices to the edge and the triangles to which the
%   intersection points belong. For this function to run properly, all
%   redundancy must be removed from the mesh (see STL_RemoveRedundancy.m).
%
%   Dependencies: TRI_IntersectWithPlane.m
%                 DistanceFromVertexToPlane.m
%                 DeleteUnreferencedElements.m
%
%   Known parents: TRI_CutWithBoundedPlane.m
%                  TRI_RegionCutter.m

%Created on 15/02/2006 by Ward Bartels.
%WB, 28/02/2006: Added calculation of missing intersection points.
%WB, 25/01/2007: Rewrite; TRI_IntersectWithPlane is now used and
%                intersecting inside a triangle is now possible.
%Stabile, fully functional.


%Calculate regular plane intersection <<TRI_IntersectWithPlane.m>>
[V_new, contour, triangles, edges_cut, distances] = TRI_IntersectWithPlane(F, V, beta);

%Calculate distances from new points to bounding planes
%<<DistanceFromVertexToPlane.m>>
dist_bound = DistanceFromVertexToPlane(V_new, beta_bound);

%Rearrange distances to bounding planes into contour edge rays
d1 = dist_bound(contour(:,1),:);
d2 = dist_bound(contour(:,2),:);

%Find rays parallel with bounding planes
parallel = abs(d2-d1)<eps;
par_inside = parallel & (d1>=0 | d2>=0); %Fail-safe when d1 and d2 differ slightly
par_outside = parallel & ~par_inside;
d1(par_inside) = Inf;                    %Will cause NaN in normcoord, ignored by max/min
d1(par_outside) = -2;
d2(par_outside) = -1;                    %Will cause 2 in normcoord, caught by max

%Check where rays run into the bounded area
entering = d1<d2;

%Calculate ray coordinates of bounding plane - ray intersections
normcoord = d1./(d1-d2);

%Find section of ray inside bounding planes
normcoord_max = normcoord;
normcoord_max(~entering) = NaN;
[max_enter, ind_enter] = max(normcoord_max, [], 2);
normcoord(entering) = NaN;
[min_exit, ind_exit] = min(normcoord, [], 2);

%Eliminate contour edges fully outside of the bounded areas
ind_delete = max_enter>=min_exit | max_enter>=1 | min_exit<=0;
contour(ind_delete,:) = []; triangles(ind_delete,:) = []; edges_cut(ind_delete,:) = [];
max_enter(ind_delete,:) = []; ind_enter(ind_delete,:) = [];
min_exit(ind_delete,:) = []; ind_exit(ind_delete,:) = [];

%Ignore contour points inside the bounding planes
max_enter(max_enter<0) = NaN;
min_exit(min_exit>1) = NaN;

%Replace triangle vertices with intersections where necessary
repl_enter = ~isnan(max_enter);
repl_exit = ~isnan(min_exit);
replace = [repl_enter repl_exit];
basept = V_new([contour(repl_enter,1); contour(repl_exit,1)],:);
edgevec = V_new([contour(repl_enter,2); contour(repl_exit,2)],:)-basept;
contour(replace) = size(V_new, 1)+(1:sum(replace(:)));
V_new = [V_new; basept+edgevec.*repmat([max_enter(repl_enter,:); min_exit(repl_exit,:)], 1, 3)];

%Remove unused vertices <<DeleteUnreferencedElements.m>>
[V_new, contour] = DeleteUnreferencedElements(V_new, contour);

%Flag replaced vertices
edges_cut(replace) = 0;

%Store "active" bounding plane indices for all contour points
failed_test = zeros(size(contour));
failed_test(replace) = [ind_enter(repl_enter); ind_exit(repl_exit)];

%Ensure modified contour points are in second column
contour(repl_enter,:) = contour(repl_enter,[2 1]);
edges_cut(repl_enter,:) = edges_cut(repl_enter,[2 1]);
failed_test(repl_enter,:) = failed_test(repl_enter,[2 1]);


% %Original code (wb20070118), begin:
% 
% %---------------------%
% % Calculate distances %
% %---------------------%
% 
% %Assemble matrix containing all edges and get triangle indices
% %<<TRI_Edges.m>>
% [edges, triangles] = TRI_Edges(F);
% 
% %Assemble column vector containing edge indices (1 for point 1 -> point 2,
% %2 for point 2 -> point 3 and 3 for point 3 -> point 1)
% edges_cut = repmat([1; 2; 3], length(triangles)/3, 1);
% 
% %Calculate distances to plane for each vertex
% %<<DistanceFromVertexToPlane.m>>
% distances = DistanceFromVertexToPlane(V, beta);
% 
% %Keep only triangles intersected at least twice
% intersected = abs(sum(sign(distances(F)), 2))<2;
% intersected = intersected(triangles);
% edges = edges(intersected,:);
% edges_cut = edges_cut(intersected);
% triangles = triangles(intersected);
% 
% %Add a filler to the end of triangles and edges_cut
% triangles(end+1) = 0;
% edges_cut(end+1) = 0;
% 
% %Stack up row indices of equal edges; the filler is an index to the last
% %element in triangles <<StackEqualElementIndices.m>>
% stack = StackEqualElementIndices(sort(edges, 2), length(triangles), 2);
% 
% %Select unique edges and reshape triangles and edges_cut so their rows
% %correspond to the rows of edges
% edges = edges(stack(:,1),:);
% edges_cut = edges_cut(stack);
% triangles = triangles(stack);
% 
% 
% %-------------------------------%
% % Calculate intersection points %
% %-------------------------------%
% 
% %Select distances from begin- and endpoints of unique edges to plane
% d1 = distances(edges(:,1));
% d2 = distances(edges(:,2));
% 
% %Limit edges, stack, d1 and d2 to edges which intersect the plane (in an
% %intersected triangle, not all edges are intersected normally); this will
% %also eliminate triangles lying in the plane
% ind = sign(d1)~=sign(d2);
% edges = edges(ind,:);
% edges_cut = edges_cut(ind,:);
% triangles = triangles(ind,:);
% d1 = d1(ind);
% d2 = d2(ind);
% 
% %Select begin- and endpoints of edges which intersect the plane
% p1 = V(edges(:,1),:);
% p2 = V(edges(:,2),:);
% 
% %Calculate intersection points for edges which intersect the plane
% V_new = (p2.*(d1*ones(1, 3))-p1.*(d2*ones(1, 3)))./((d1-d2)*ones(1, 3));
% 
% %Check which intersection points fail which boundary test
% %<<DistanceFromVertexToPlane.m>>
% distances_bound = DistanceFromVertexToPlane(V_new, beta_bound);
% outside_bound = [true(size(distances_bound, 1), 1) distances_bound<0] ;
% [ignoble, failed_test] = max(cumsum(outside_bound, 2), [], 2);
% failed_test = failed_test-1;
% 
% %Assemble vector of indices into V_new; rearrange triangles, failed_test
% %and edges_cut
% contour = [(1:size(triangles, 1))'; find(triangles(:,2)~=0)];
% failed_test = failed_test(contour);
% triangles = nonzeros(triangles);
% edges_cut = nonzeros(edges_cut);
% 
% 
% %------------------------------------%
% % Remove invalid intersection points %
% %------------------------------------%
% 
% %Remove redundancy from V_new
% [ignoble, ind, jnd] = unique(single(V_new), 'rows');  %Single precision for thresholding
% V_new = V_new(ind,:);
% contour = jnd(contour);
% 
% %Remove duplicate points occurring when a triangle is cut through a corner
% %point
% [triangles_contour, ind] = unique([triangles contour], 'rows');
% triangles = triangles_contour(:,1);
% contour   = triangles_contour(:,2).';
% edges_cut   = edges_cut(ind).';
% failed_test = failed_test(ind).';
% 
% %Stack up equal triangles <<StackEqualElementIndices.m>>
% stack = StackEqualElementIndices(triangles, 0, 2);
% 
% %Eliminate triangles cut only once
% stack(stack(:,2)==0,:) = [];
% 
% %Assemble matrices representing contour edges, cut edges and failed
% %boundary tests; select unique triangles
% triangles   = triangles(stack(:,1));
% contour     = contour(stack);
% edges_cut   = edges_cut(stack);
% failed_test = failed_test(stack);
% 
% %Remove triangles for which no intersection points pass the boundary test
% failed = failed_test~=0;
% failed_all = all(failed, 2);
% contour(failed_all,:)     = [];
% edges_cut(failed_all,:)   = [];
% triangles(failed_all,:)   = [];
% failed_test(failed_all,:) = [];
% failed(failed_all,:)      = [];
% 
% %Remove all intersection points which fail a boundary test
% contour(failed)   = 0;
% edges_cut(failed) = 0;
% [failed_test, ind] = max(failed_test, [], 2);
% ind_flip = ind==1;
% contour(ind_flip,:)   = contour(ind_flip,[2 1]);
% edges_cut(ind_flip,:) = edges_cut(ind_flip,[2 1]);
% 
% %Remove unused vertices <<DeleteUnreferencedElements.m>>
% [V_new, contour] = DeleteUnreferencedElements(V_new, contour, 0);
% 
% 
% %---------------------------------------%
% % Calculate missing intersection points %
% %---------------------------------------%
% 
% %Find triangles cut in only one edge by half-plane
% ind_half = ~contour(:,2);
% 
% %Build intersection plane matrices
% ind = failed_test(ind_half)*3;
% beta1 = beta_bound([ind-2 ind-1 ind]',:);
% beta2 = repmat(beta, length(ind), 1);
% beta3 = V(F(triangles(ind_half),:)',:);
% 
% %Calculate intersection points <<IntersectThreePlanes.m>>
% V_inter = IntersectThreePlanes(beta1, beta2, beta3);
% 
% %Reference intersection points and add them to V_new
% contour(ind_half,2) = (size(V_new, 1)+1:size(V_new, 1)+size(V_inter, 1)).';
% V_new = [V_new; V_inter];
% 
% %Original code (wb20070118), end.