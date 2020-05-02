function conn = Graph_Edges2Connectivity(edges, form, directed, loops, numvert)
%Graph_Edges2Connectivity  Generate connectivity array from edge list.
%
%   Syntax:
%    conn = Graph_Edges2Connectivity(edges, form, directed, loops, numvert)
%
%   Input:
%    edges:    N-by-2 array containing edges. Each row represents an edge,
%              each element is a link to a vertex. If the graph is
%              directed, edges will run from the first column to the
%              second. If it is not directed, only one edge may be provided
%              per connected vertex pair.
%    form:     String indicating which form the output connectivity array
%              should take. Can be set to 'logsparse' (square logical
%              sparse adjacency matrix), 'cell' (cell array adjacency list)
%              or 'sparse' (sparse adjacency list).
%    directed: Logical indicating whether or not the graph is directed.
%              Optional, defaults to false.
%    loops:    Logical indicating whether or not loop edges should be added
%              to all vertices. Optional, defaults to false.
%    numvert:  Total number of vertices in the graph. Optional, defaults to
%              max(edges(:)).
%
%   Output:
%    conn: Connectivity array with shape depending on the string passed to
%          form. If form is set to...
%           'logsparse', it is an N-by-N logical sparse adjacency matrix.
%               Each row and column correspond to a vertex and each true
%               element indicates an edge from its row to its column.
%           'cell', it is an N-by-1 cell array adjacency list containing
%               column vectors of vertex indices. Each column vector
%               corresponds with a vertex and lists all vertices connected
%               to it.
%           'sparse', it is an M-by-N sparse array adjacency list
%               containing vertex indices. Each row corresponds with a
%               vertex and lists all vertices connected to it.
%
%   Effect: This function will return a connectivity matrix of the
%   requested shape, based on the provided edge list which describes a
%   graph. It is the opposite of Graph_Connectivity2Edges.m.
%
%   Dependencies: IncrementalRuns.m
%
%   Known parents: ICP_MatchDatasets.m
%                  Contour_Loops.m
%                  TRI_VertexNeighbourhood.m
%                  Graph_Components.m
%                  Graph_BreadthFirstTree.m
%                  SkeletalModel/private/Update.m

%Created on 22/12/2006 by Ward Bartels.
%WB, 05/07/2007: Added directed and loops arguments.
%WB, 13/07/2007: Speed improvements for logsparse.
%Stabile, fully functional.


%Set defaults for directed and loops
directed = nargin>=3 && directed;
loops = nargin>=4 && loops;

%Calculate number of vertices if not provided
if nargin<5
    numvert = max(max(edges));
end

%Check if logical sparse was requested
if strcmpi(form, 'logsparse')
    
    %Build square logical sparse connectivity matrix
    conn = sparse(edges(:,1), edges(:,2), true, numvert, numvert);
    
    %Add inverted edges if the graph is undirected
    if ~directed
        conn = conn | conn.';
    end
    
    %Add loop edges if requested
    if loops
        conn = conn | speye(size(conn, 1));
    end
    
else
    
    %Add inverted edges if the graph is undirected
    if ~directed
        edges = [edges; edges(:,[2 1])];
    end
    
    %Add loop edges if requested
    if loops
        edges = [edges; (1:numvert).'*[1 1]];
    end
    
    %Sort edges
    edges = sortrows(edges);
    
    %Check if sparse or cell array was requested
    if strcmpi(form, 'cell')
        
        %Build connectivity cell array
        conn = accumarray(edges(:,1), edges(:,2), [numvert 1], @(x) {x});
        
    elseif strcmpi(form, 'sparse')
        
        %Build sparse connectivity matrix <<IncrementalRuns.m>>
        cols = IncrementalRuns(accumarray(edges(:,1), 1));
        conn = sparse(edges(:,1), cols, edges(:,2), numvert, max(cols));
        
    else
        error([mfilename ':UnknownForm'], 'Unrecognised form string.');
    end
end