function [V_new, contour, triangles, edges_cut, distances] = TRI_IntersectWithPlane(F, V, beta)
%TRI_IntersectWithPlane  Intersect mesh and plane.
%
%   Syntax:
%    [V_new, contour, triangles, edges_cut, distances] = ...
%        TRI_IntersectWithPlane(F, V, beta)
%
%   Input:
%    F:    N-by-3 array containing indices into V. Each row represents a
%          triangle, each element is a link to a vertex in V.
%    V:    N-by-3 array containing vertex coordinates. Each row represents
%          a vertex; the first, second and third columns represent X-, Y-
%          and Z-coordinates respectively.
%    beta: 3-by-3 array which defines the plane by specifying vertex
%          coordinates. Each row represents a vertex that lies in the
%          plane; the first, second and third colums represent X-, Y- and
%          Z-coordinates respectively.
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
%                 columns 2 and 3, and so on.
%    distances:   Column vector containing distances from each vertex to
%                 the specified plane. Each row corresponds to the same row
%                 in V.
%
%   Effect: This function will extract all edges from the mesh and
%   determine which edges intersect the plane. For these edges, the
%   intersection points are calculated. The coordinates of these points are
%   returned along with indices to the edge and the triangles to which the
%   intersection points belong. For this function to run properly, all
%   redundancy must be removed from the mesh (see STL_RemoveRedundancy.m).
%
%   Dependencies: TRI_Edges.m
%                 DistanceFromVertexToPlane.m
%                 StackEqualElementIndices.m
%
%   Known parents: TRI_CutWithPlane.m
%                  TRI_IntersectWithBoundedPlane.m

%Created on 19/12/2005 by Ward Bartels.
%WB, 03/01/2006: Added removal of triangles intersected in a single point;
%                added additional outputs.
%WB, 04/01/2006: Removed output of triangles; performance improvements.
%WB, 06/01/2006: Performance improvements.
%WB, 03/01/2007: Rewrite; outputs are now similar to those in
%                TRI_IntersectWithBoundedPlane.m.
%WB, 30/04/2007: Performance improvements.
%WB, 24/08/2007: Rewrite; better robustness and slight speed increase.
%Stabile, fully functional.


%Get mesh edges <<TRI_Edges.m>>
edges = TRI_Edges(F);

%Calculate distance to plane for vertices <<DistanceFromVertexToPlane.m>>
distances = DistanceFromVertexToPlane(V, beta);

%Determine which vertices are above the plane (positive normal)
above = distances>0; %plane is virtually "moved upwards" by >

%Determine which edges are intersected
cut = xor(above(edges(:,1)), above(edges(:,2)));

%Return if no triangles are cut
if ~any(cut)
    V_new = zeros(0, 3);
    contour = zeros(0, 2);
    triangles = zeros(0, 1);
    edges_cut = zeros(0, 2);
    return
end

%Keep only intersected edges
edges = edges(cut,:);

%Isolate unique edges <<StackEqualElementIndices.m>>
ind = [find(cut); 0];
stack = StackEqualElementIndices(sort(edges, 2), numel(ind));
edges = edges(stack(:,1),:);
stack = ind(stack); %the original edge indices in stack can tell which triangles were intersected
%the above could also be done before finding intersected edges; however,
%it's usually faster to first remove non-intersecting edges and then look
%for duplicate edges (which is slower and thus benefits from having fewer
%edges to deal with)

%Isolate data of intersected edges (vertices and distance to plane)
p1 = V(edges(:,1),:);
p2 = V(edges(:,2),:);
d1 = distances(edges(:,1));
d2 = distances(edges(:,2));

%Calculate intersection points for edges that intersect the plane
dtot = d1-d2;
V_new = p2.*((d1./dtot)*ones(1, 3))-p1.*((d2./dtot)*ones(1, 3));

%Determine intersected triangles and triangle sides
[contour, ignoble, stack] = find(stack);
edges_cut = mod(stack-1, 3)+1;
triangles = (stack-edges_cut)/3+1;

%Consolidate intersected triangles <<StackEqualElementIndices.m>>
stack = StackEqualElementIndices(triangles);
contour = contour.';
contour = contour(stack);
edges_cut = edges_cut.';
edges_cut = edges_cut(stack);
triangles = triangles(stack(:,1));


% %Original code (wb20070824), begin:
% 
% %---------------------%
% % Calculate distances %
% %---------------------%
% 
% %Assemble matrix containing all edges <<TRI_Edges.m>>
% edges = TRI_Edges(F);
% 
% %Calculate distances to plane for each vertex
% %<<DistanceFromVertexToPlane.m>>
% distances = DistanceFromVertexToPlane(V, beta);
% 
% %List triangles intersected at least twice, and their edges
% intersected = abs(sum(sign(distances(F)), 2))<2;
% triangles = find(intersected);
% edges_cut = [1;2;3]*ones(1, numel(triangles));
% edges_cut = edges_cut(:);
% triangles = triangles.';
% triangles = triangles([1 1 1],:);
% triangles = triangles(:);
% intersected = intersected.';
% edges = edges(intersected([1 1 1],:),:);
% 
% %Add a filler to the end of triangles and edges_cut
% triangles(end+1) = 0;
% edges_cut(end+1) = 0;
% 
% %Stack up row indices of equal edges; the filler is an index to the last
% %element in triangles <<StackEqualElementIndices.m>>
% stack = StackEqualElementIndices(sort(edges, 2), numel(triangles), 2);
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
% %Assemble vector of indices into V_new; rearrange triangles and edges_cut
% contour = [(1:size(triangles, 1))'; find(triangles(:,2)~=0)];
% triangles = nonzeros(triangles);
% edges_cut = nonzeros(edges_cut);
% 
% 
% %------------------------------------%
% % Remove invalid intersection points %
% %------------------------------------%
% 
% %Remove redundancy from V_new
% [ignoble, ind, jnd] = unique(V_new, 'rows'); %single(V_new) removed, experimental!
% V_new = V_new(ind,:);
% contour = jnd(contour);
% 
% %Remove duplicate points occurring when a triangle is cut through a corner
% %point
% [triangles_contour, ind] = unique([triangles contour], 'rows');
% triangles = triangles_contour(:,1);
% contour   = triangles_contour(:,2).';
% edges_cut = edges_cut(ind).';
% 
% %Stack up equal triangles <<StackEqualElementIndices.m>>
% stack = StackEqualElementIndices(triangles, 0, 2);
% 
% %Eliminate triangles cut only once
% stack(stack(:,2)==0,:) = [];
% 
% %Assemble matrices representing contour edges and cut edges; select unique
% %triangles
% triangles = triangles(stack(:,1));
% contour   = contour(stack);
% edges_cut = edges_cut(stack);
% 
% %Remove unused vertices <<DeleteUnreferencedElements.m>>
% [V_new, contour] = DeleteUnreferencedElements(V_new, contour, 0);
% 
% %Original code (wb20070824), end.


% %Original code (wb20070103), begin:
% 
% %---------------------%
% % Calculate distances %
% %---------------------%
% 
% %Assemble matrix containing all edges and get triangle indices
% %<<TRI_Edges.m>>
% [edges, triangles] = TRI_Edges(F);
% 
% %Calculate distances to plane for each vertex
% %<<DistanceFromVertexToPlane.m>>
% distances = DistanceFromVertexToPlane(V, beta);
% 
% %Keep only triangles intersected at least twice
% intersected = abs(sum(sign(distances(F)), 2))<2;
% edges = edges(intersected(triangles),:);
% triangles = triangles(intersected(triangles));
% 
% %Add a filler to the end of triangles
% triangles(end+1) = 0;
% 
% %Stack up row indices of equal edges; the filler is an index to the last
% %element in triangles <<StackEqualElementIndices.m>>
% stack = StackEqualElementIndices(sort(edges, 2), length(triangles), 2);
% 
% %Select unique edges and reshape triangles so its rows correspond to the
% %rows of edges
% edges = edges(stack(:,1),:);
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
% 
% %------------------------------------%
% % Remove invalid intersection points %
% %------------------------------------%
% 
% %Assemble vector of indices into edges and V_new; rearrange triangles
% pointer_V_new = [[1:size(edges, 1)]'; find(triangles(:,2)~=0)];
% pointer_edges = pointer_V_new;
% triangles = nonzeros(triangles);
% 
% %Remove redundancy from V_new
% [ignoble, ind, jnd] = unique(single(V_new), 'rows'); %Single precision for thresholding
% V_new = V_new(ind,:);
% pointer_V_new = jnd(pointer_V_new);
% 
% %Remove duplicate points occurring when a triangle is cut through a corner
% %point <<DeleteUnreferencedElements.m>>
% [triangles_pointer_V_new, ind] = unique([triangles pointer_V_new], 'rows');
% triangles = triangles_pointer_V_new(:,1);
% pointer_V_new = triangles_pointer_V_new(:,2).';
% [V_new, pointer_V_new] = DeleteUnreferencedElements(V_new, pointer_V_new);
% 
% %Stack up equal triangles <<StackEqualElementIndices.m>>
% stack = StackEqualElementIndices(triangles, 0, 2);
% stack(~stack(:,2),:) = []; %Fail-safe
% 
% %Assemble matrix representing contour edges
% edges_contour = pointer_V_new(stack);
% 
% %Distinguish between output of contour only or output of full data
% if nargout<=2 %Output of contour only (e.g. reslicing)
%     
%     %Remove redundant contour edges occurring when a triangle is cut
%     %through a side (two triangles will be cut, and will both contribute
%     %the same edge to the contour)
%     edges_contour = unique(edges_contour, 'rows');
%     
% else %Output of full data (e.g. cutting a mesh)
%     
%     %Remove edge pointers for duplicate points
%     %<<DeleteUnreferencedElements.m>>
%     pointer_edges = pointer_edges(ind);
%     [edges, pointer_edges] = DeleteUnreferencedElements(edges, pointer_edges);
%     
%     %Assemble matrices representing intersected edges
%     edge1 = edges(pointer_edges(stack(:,1)),:);
%     edge2 = edges(pointer_edges(stack(:,2)),:);
% end
% 
% %Original code (wb20070103), end.