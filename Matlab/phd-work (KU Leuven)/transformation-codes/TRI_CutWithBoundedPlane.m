function [F, V, contour] = TRI_CutWithBoundedPlane(varargin)
%TRI_CutWithBoundedPlane  Cut mesh by intersecting with a bounded plane.
%wb20070823
%
%   Syntax:
%    [F, V, contour] = TRI_CutWithBoundedPlane(F, V, beta, beta_bound)
%    [F, V, contour] = TRI_CutWithBoundedPlane(F, V, V_new, contour, ...
%                                              triangles, edges_cut)
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
%    V_new:      N-by-3 array containing vertex coordinates of inter-
%                section points. Each row represents a vertex; the first,
%                second and third columns represent X-, Y- and Z-
%                coordinates respectively.
%    contour:    N-by-2 array containing the intersection line. Each row
%                represents an edge and corresponds to the same row in
%                triangles; each element is a row index into V_new.
%    triangles:  Column vector containing the intersected triangles. Each
%                element is a row index into F.
%    edges_cut:  N-by-2 array indicating which edge of a triangle is
%                intersected by the plane. Each row corresponds to the same
%                row in triangles and each column corresponds to the same
%                column in contour. A 1 indicates the edge connecting
%                columns 1 and 2 in F, a 2 indicates the edge connecting
%                columns 2 and 3, and so on. For triangles intersected in
%                only one edge, edges_cut should contain a zero in the
%                second column; both columns should contain a zero when
%                intersecting inside a triangle.
%
%   Output:
%    F:        N-by-3 array containing indices into V. Each row represents
%              a triangle, each element is a link to a vertex in V.
%    V:        N-by-3 array containing vertex coordinates. Each row
%              represents a vertex; the first, second and third columns
%              represent X-, Y- and Z-coordinates respectively.
%    contour:  N-by-2 array containing the intersection line. Each row
%              represents an edge, each element is a row index into V.
%
%   Effect: This function will determine the intersection line of a mesh
%   and a bounded plane. The intersection line is added to the mesh, and
%   edges lying along the intersection are returned. This function may
%   introduce unused vertices into V. The second syntax can be used to plug
%   the (possibly modified) output of a call to
%   TRI_IntersectWithBoundedPlane.m directly into this function.
%
%   Dependencies: TRI_IntersectWithBoundedPlane.m
%                 TRI_Normals.m
%                 TRI_Edges.m
%                 StackEqualElementIndices.m
%                 Graph_Components.m
%
%   Known parents: TRI_CutWithMultiPlane.m
%                  TRI_RegionCutter.m

%Created on 01/03/2006 by Ward Bartels.
%WB, 02/03/2006: Added parsing of output variables, removed duplication of
%                intersection line.
%WB, 04/01/2007: Simplified triangle assembling code.
%WB, 17/01/2007: Cutting inside a single triangle is now possible.
%WB, 30/04/2007: Speed improvements.
%WB, 23/08/2007: Added handling of zero-length edges.
%Stabile, fully functional.


%--------------------------%
% Intersect plane and mesh %
%--------------------------%

%Intersect mesh with plane if necessary <<TRI_IntersectWithBoundedPlane.m>>
if nargin==4
    [F, V, beta, beta_bound] = varargin{:};
    [V_new, contour, triangles, edges_cut] = TRI_IntersectWithBoundedPlane(F, V, beta, beta_bound);
else
    [F, V, V_new, contour, triangles, edges_cut] = varargin{:};
end

%Add V_new to V and reference it
contour = contour+size(V, 1);
V = [V; V_new];

%Split off triangles cut on the inside (not through the edges)
ind_i = edges_cut(:,1)==0;
contour_i   = contour(ind_i,:);   contour(ind_i,:)   = [];
triangles_i = triangles(ind_i,:); triangles(ind_i,:) = [];
edges_cut(ind_i,:) = [];

%Distinguisch between fully cut and half cut triangles
ind_h = edges_cut(:,2)==0;
contour_h   = contour(ind_h,:);   contour(ind_h,:)   = [];
triangles_h = triangles(ind_h,:); triangles(ind_h,:) = [];
edges_cut_h = edges_cut(ind_h,1); edges_cut(ind_h,:) = [];

%Make sure edges_cut is in "ascending" order horizontally (1-2, 2-3, 3-1)
ind_flip = logical(mod(diff(edges_cut, 1, 2)+2, 3));
edges_cut(ind_flip,:) = edges_cut(ind_flip,[2 1]);
contour(ind_flip,:)   = contour(ind_flip,[2 1]);


%---------------------------------------------------%
% Subdivide triangles intersected through two edges %
%---------------------------------------------------%

%Find common ("tip") and uncommon ("base") points between intersected edges
tip   = F(sub2ind(size(F), triangles, edges_cut(:,2)));
base1 = F(sub2ind(size(F), triangles, edges_cut(:,1)));
base2 = F(sub2ind(size(F), triangles, 6-sum(edges_cut, 2)));

%Fail-safe
if isempty(tip), tip = zeros(0, 1); end
if isempty(base1), base1 = zeros(0, 1); end
if isempty(base2), base2 = zeros(0, 1); end

%Assemble new F
F_new1 = [base1(:,1)   contour(:,1) base2(:,1);
          contour(:,1) contour(:,2) base2(:,1);
          contour(:,1) tip(:,1)     contour(:,2)];


%--------------------------------------------------%
% Subdivide triangles intersected through one edge %
%--------------------------------------------------%

%Find points adjacent and non-adjacent to the intersected edge
nonadjacent = F(sub2ind(size(F), triangles_h, mod(edges_cut_h+1, 3)+1));
adjacent1   = F(sub2ind(size(F), triangles_h, mod(edges_cut_h, 3)+1));
adjacent2   = F(sub2ind(size(F), triangles_h, edges_cut_h));

%Fail-safe
if isempty(nonadjacent), nonadjacent = zeros(0, 1); end
if isempty(adjacent1), adjacent1 = zeros(0, 1); end
if isempty(adjacent2), adjacent2 = zeros(0, 1); end

%Assemble new F
F_new2 = [contour_h(:,1)   adjacent1(:,1) contour_h(:,2);
          contour_h(:,1)   contour_h(:,2) adjacent2(:,1);
          contour_h(:,2)   adjacent1(:,1) nonadjacent(:,1);
          nonadjacent(:,1) adjacent2(:,1) contour_h(:,2)];


%-----------------------------------------------%
% Subdivide triangles intersected on the inside %
%-----------------------------------------------%

%Create vectors pointing from first contour point to triangle corner points
corners = V(F(triangles_i,:),:)-V(repmat(contour_i(:,1), 3, 1),:);

%Check between which two corner vectors the second contour point lies
%<<TRI_Normals.m>>
cp = cross(corners, repmat(V(contour_i(:,2),:)-V(contour_i(:,1),:), 3, 1), 2);
cp = reshape(dot(cp, repmat(TRI_Normals(F(triangles_i,:), V, false), 3, 1), 2)>0, [], 3);
distant = ones(size(triangles_i)); %Between vectors 2 and 3
distant(cp(:,3) & ~cp(:,1)) = 2;   %Between vectors 3 and 1
distant(cp(:,1) & ~cp(:,2)) = 3;   %Between vectors 1 and 2

%Distribute points based on corner vectors
corner1 = F(sub2ind(size(F), triangles_i, distant));
corner2 = F(sub2ind(size(F), triangles_i, mod(distant, 3)+1));
corner3 = F(sub2ind(size(F), triangles_i, mod(distant+1, 3)+1));

%Fail-safe, in case F has only one row
corner1 = corner1(:);
corner2 = corner2(:);
corner3 = corner3(:);

%Assemble new F
F_new3 = [corner3 corner1        contour_i(:,1);
          corner1 corner2        contour_i(:,1);
          corner2 corner3        contour_i(:,2);
          corner3 contour_i(:,1) contour_i(:,2);
          corner2 contour_i(:,2) contour_i(:,1)];


%--------------------------%
% Remove zero-length edges %
%--------------------------%

%Remove triangles intersected by the plane
F([triangles; triangles_h; triangles_i],:) = [];

%Compile contour lines
contour = [contour; contour_h; contour_i];

%Assemble all new triangles
F_new = [F_new1; F_new2; F_new3];

%Find unique edges in new triangles <<TRI_Edges.m>>
%                                   <<StackEqualElementIndices.m>>
edges = TRI_Edges(F_new);
stack = StackEqualElementIndices(sort(edges, 2));
edges = edges(stack(:,1),:);

%Check which edges are zero length
degenerate = all(single(V(edges(:,1),:))==single(V(edges(:,2),:)), 2);

%If any edges are degenerate, remove triangles
if any(degenerate)
    
    %Remove new triangles with faulty edges
    F_new(ceil(nonzeros(stack(degenerate,:))/3),:) = [];
    
    %Check which groups of edge points coincide <<Graph_Components.m>>
    edges = edges(degenerate,:);
    [points, ignoble, jnd] = unique(edges);
    group = Graph_Components(reshape(jnd, [], 2)); %group number for each point
    
    %Determine which points must be replaced by which
    [group, ind] = sort(group);
    points = points(ind);
    first = [true; group(1:end-1)~=group(2:end)]; %true for first group element
    newpoints = points(first);
    newpoints = newpoints(group(~first));
    points = points(~first); %now replace points by corresponding newpoints
    
    %Weld coinciding vertices in new triangles
    [repl, loc] = ismember(F_new, points);
    F_new(repl) = newpoints(loc(repl));
    [repl, loc] = ismember(contour, points);
    contour(repl) = newpoints(loc(repl));
    
    %Weld coinciding vertices in old triangles (should be rare, since
    %points is preferentially high and thus inside V_new)
    ind = points>size(V, 1)-size(V_new, 1);
    if ~all(ind)
        points(ind,:) = [];
        newpoints(ind,:) = [];
        [repl, loc] = ismember(F, points);
        F(repl) = newpoints(loc(repl));
    end
end

%Add new triangles
F = [F; F_new];