function components = Graph_Components(argin, varargin)
%Graph_Components  Isolate connected components of an undirected graph.
%
%   Syntax:
%    components = Graph_Components(conn)
%    components = Graph_Components(edges, numvert)
%
%   Input:
%    conn:    N-by-N logical sparse adjacency matrix. Each row and column
%             correspond to a vertex and each true element indicates an
%             edge from its row to its column.
%    edges:   N-by-2 array containing edges. Each row represents an edge,
%             each element is a link to a vertex. Only one edge may be
%             provided per connected vertex pair.
%    numvert: Total number of vertices in the graph. Optional, defaults to
%             max(edges(:)).
%
%   Output:
%    components: Column vectors containing component numbers. Each element
%                corresponds to a vertex and indicates to which component
%                the corresponding vertex belongs.
%
%   Effect: This function will isolate connected components from a graph,
%   represented by its adjacency matrix. It assumes the provided graph is
%   undirected.
%
%   Dependencies: Graph_Edges2Connectivity.m
%                 RunLengthDecode.m
%
%   Known parents: TRI_SeparateShells.m
%                  TRI_CutWithBoundedPlane.m
%                  SkeletalModel/private/Update.m

%Created on 12/07/2007 by Ward Bartels.
%Stabile, fully functional.


%Ensure correct type of input argument <<Graph_Edges2Connectivity.m>>
if islogical(argin) %first syntax, add loop edges to adjacency matrix
    conn = argin | speye(size(argin, 1));
else                %second syntax, build adjacency matrix from edge list
    conn = Graph_Edges2Connectivity(argin, 'logsparse', false, true, varargin{:});
end

%Find connected components
[perm, ignoble, blocks] = sparsfun('dmperm', conn); %inlined dmperm.m

%Give each vertex a component label <<RunLengthDecode.m>>
components = RunLengthDecode(diff(blocks.'));
components(perm,1) = components;