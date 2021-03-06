function edges = Graph_Simplices2Edges(simplices)
%Graph_Simplices2Edges  Generate edge list from simplices.
%
%   Syntax:
%    edges = Graph_Simplices2Edges(simplices)
%
%   Input:
%    simplices: M-by-N array containing tessellation simplices (triangles,
%               tetrahedrons or (N-1)-dimensional hypertriangles). Each row
%               represents a simplex, each element is a link to a vertex.
%
%   Output:
%    edges: L-by-2 array containing edges. Each row represents an edge,
%           each element is a link to a vertex. Each edge in this array
%           will occur exactly once.
%
%   Effect: From a list of simplices, for example generated by delaunay.m
%   or delaunayn.m, this function will generate an edge list that describes
%   the corresponding graph. In case of a triangulation's face matrix, it
%   is more efficient to use TRI_Edges.m.
%
%   Dependencies: none
%
%   Known parents: ICP_MatchDatasets.m

%Created on 19/12/2006 by Ward Bartels.
%Stabile, fully functional.


%Determine all possible combinations of vertices in one simplex
ind = combnk(1:size(simplices, 2), 2);

%List all edges
edges = reshape(simplices(:,ind(:)), [], 2);

%Remove redundant edges
edges = unique(sort(edges, 2), 'rows');