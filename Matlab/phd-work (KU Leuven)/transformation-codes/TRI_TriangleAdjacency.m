function [connections, side, edges] = TRI_TriangleAdjacency(argin, filler)
%TRI_TriangleAdjacency  Determine which triangles share an edge.
%
%   Syntax:
%    [connections, side, edges] = TRI_TriangleAdjacency(F, filler)
%    [connections, side, edges] = TRI_TriangleAdjacency(edges, filler)
%
%   Input:
%    F:      M-by-3 array defining a triangle mesh. The rows correspond to
%            different triangles and the columns correspond to the three
%            vertices that make up each triangle. The elements are vertex
%            indices.
%    edges:  First output of TRI_Edges(F); rows can be set to NaN to turn
%            corresponding edges into border edges.
%    filler: Scalar that will be used to fill the empty spaces in
%            connections. Optional, defaults to 0.
%
%   Output:
%    connections: N-by-2 array indicating which triangles are connected.
%                 The rows correspond to unique edges and the columns
%                 correspond to the two triangles that are connected by
%                 each edge. The elements are row indices into F. For
%                 border edges, the specified filler will appear in the
%                 second column.
%    side:        N-by-2 array indicating on which side the triangles are
%                 connected. The rows correspond to unique edges and the
%                 columns correspond to the two triangles that are
%                 connected by each edge. The elements can be 1, 2 and 3,
%                 with a 1 indicating that the corresponding triangle is
%                 connected on its first edge (point 1 and 2), and so on.
%                 For border edges, the same fillers are used as in
%                 connections.
%    edges:       N-by-2 array containing connecting edges. The rows
%                 correspond to unique edges and the columns correspond to
%                 the two vertices that make up each edge. The elements are
%                 vertex indices.
%
%   Effect: This function will produce a list of which triangles in the
%   provided mesh are connected by a shared edge. If more than two
%   triangles have an edge in common, a warning will be produced, and all
%   but two triangles will be listed as containing a border edge rather
%   than sharing the same edge.
%
%   Dependencies: TRI_Edges.m
%                 StackEqualElementIndices.m
%
%   Known parents: TRI_SeparateShells.m
%                  TRI_DetermineTriangleConnections.m
%                  TRI_RemoveInvertedNormals.m
%                  TRI_Isthmus.m

%Created on 03/01/2008 by Ward Bartels.
%WB, 16/07/2008: Added output of edges.
%Stabile, fully functional.


%Extract edges if necessary <<TRI_Edges.m>>
if size(argin, 2)>2
    edges = TRI_Edges(argin);
else
    edges = argin;
end

%Stack up row indices of equal edges <<StackEqualElementIndices.m>>
stack = StackEqualElementIndices(sort(edges, 2));

%Handle mesh containing only border edges
if size(stack, 2)==1
    stack(:,2) = 0;
end

%Handle faulty meshes
faulty = nonzeros(stack(:,3:end));
if ~isempty(faulty)
    warning([mfilename ':FaultyMesh'], 'More than 2 triangles share an edge');
    stack = [stack(:,[1 2]); faulty zeros(size(faulty))];
end

%Translate edge indices into triangle indices (relies on TRI_Edges.m)
connections = ceil(stack/3);

%Get indices (1-3) of connected edge
if nargout>=2
    side = mod(stack-1, 3)+1;
    side(connections==0) = 0;
end

%Get connecting (unique) edges
if nargout>=3
    edges = edges(stack(:,1),:);
end

%Change filler if necessary
if nargin>=2
    connections(connections==0) = filler;
    if nargout>=2
        side(side==0) = filler;
    end
end