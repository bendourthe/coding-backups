function [from, to] = Contour_VertexConnections(argin)
%Contour_VertexConnections  Determine which vertices are connected.
%
%   Syntax:
%    [from, to] = Contour_VertexConnections(edges)
%    [from, to] = Contour_VertexConnections(F)
%
%   Input:
%    edges: N-by-2 array containing edges. Each row represents an edge,
%           each element is a link to a vertex. Each edge in this array
%           must occur exactly once.
%    F:     N-by-3 array containing vertex indices. Each row represents a
%           triangle, each element is a link to a vertex.
%
%   Output:
%    from: Column vector containing vertex indices. Each element is a link
%          to a vertex from which one or more edges "depart". This variable
%          is sorted in ascending order.
%    to:   M-by-N array containing vertex indices. Each row contains links
%          to all the vertices that share an edge with (and, therefore, are
%          connected to) the vertex on the same row in from. This array
%          list the "arrival" vertices for the edges departing in from.
%          This array will be padded with zeros if necessary; if the second
%          syntax is used, it will be sparse.
%
%   Effect: This function will find connections between vertices by looking
%   for edges which share common vertices. This function will only run
%   properly if redundancy has been removed from edges or F.
%
%   Dependencies: TRI_Edges.m
%                 StackEqualElementIndices.m
%
%   Known parents: TRI_RegionCutter.m

%Created on 06/01/2006 by Ward Bartels.
%Stabile, fully functional.


%Process input
sparse_out = size(argin, 2)>2;
if sparse_out %Second syntax: input is face array <<TRI_Edges.m>>
    edges = TRI_Edges(argin);
else          %First syntax: input is edge array
    edges = [argin; fliplr(argin)];
end

%Stack up equal points; pad empty elements with a ref to a filler element
%to be added later <<StackEqualElementIndices.m>>
stack = StackEqualElementIndices(edges(:,1), size(edges, 1)+1, Inf, sparse_out);

%Assemble from- and to-arrays
from = edges(stack(:,1),1);
if sparse_out
    [ind, jnd, s] = find(stack); %Find is very fast for sparse matrices
    to = sparse(ind, jnd, edges(s,2));
else
    to = [edges(:,2); 0];
    to = to(stack);
    if ~isequal(size(to), size(stack)), to = to.'; end %Ensure to has the correct shape
end