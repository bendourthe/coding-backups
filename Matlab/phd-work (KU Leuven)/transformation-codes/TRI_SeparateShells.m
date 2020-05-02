function [shells, F_new, varargout] = TRI_SeparateShells(varargin)
%TRI_SeparateShells  Separates a mesh into shells.
%
%   Syntax:
%    [shells, F_new, V_new, closed] = TRI_SeparateShells(F, V, edges)
%    [shells, F_new, closed] = TRI_SeparateShells(F, edges)
%    shells = TRI_SeparateShells(edges)
%
%   Input:
%    F:      N-by-3 array containing indices into V. Each row represents a
%            triangle, each element is a link to a vertex in V.
%    V:      N-by-3 array containing vertex coordinates. Each row
%            represents a  vertex; the first, second and third columns
%            represent X-, Y- and Z-coordinates respectively.
%    edges:  First output of TRI_Edges.m; rows can be set to NaN to turn
%            corresponding edges into border edges. Optional if F is
%            provided.
%
%   Output:
%    shells: Column vector containing shell numbers. Each element
%            corresponds to a triangle and a row in F. The numbers in this
%            vector indicate to which shell each triangle belongs.
%    F_new:  N-by-1 cell array containing separate F-matrices for each
%            shell. If the first syntax is used, F_new will reference into
%            V_new; otherwise, F_new will reference into V.
%    V_new:  N-by-1 cell array containing separate V-matrices for each
%            shell.
%    closed: Column vector of logicals indicating whether or not the
%            corresponding shell is a closed mesh.
%
%   Effect: This function will split up a mesh into separate shells. Shells
%   consist of a set of triangles interconnected by edges (not points).
%
%   Dependencies: TRI_Edges.m
%                 TRI_TriangleAdjacency.m
%                 Graph_Components.m
%                 DeleteUnreferencedElements.m
%
%   Known parents: TRI_SplitWithMultiPlane.m
%                  TRI_MeshCutter.m
%                  TRI_CutWithDisplacedContour.m
%                  Muscle_CutByRegionBorder.m
%                  Muscle_SelectByRegionBorder.m
%                  TRI_FlipNormalsToConvex.m

%Created on 24/01/2006 by Ward Bartels.
%WB, 15/03/2006: Added second syntax.
%WB, 10/05/2006: Replaced two for loops with cellfun calls.
%WB, 05/12/2006: Rewrite based on MatlabBGL; removed third syntax.
%WB, 13/07/2007: Replaced MatlabBGL with Graph_Components.
%WB, 03/01/2007: Moved code to TRI_TriangleAdacency.m and added syntax.
%Stabile, fully functional.


%Handle input
edges_in = size(varargin{nargin}, 2)==2;
if edges_in
    edges = varargin{nargin};
    varargin(:,end) = [];
end
if numel(varargin)>=1, F = varargin{1}; end
V_in = numel(varargin)>=2;
if V_in, V = varargin{2}; end

%Determine edges if necessary <<TRI_Edges.m>>
if ~edges_in
    edges = TRI_Edges(F);
end

%Determine which triangles share an edge <<TRI_TriangleAdjacency.m>>
conn = TRI_TriangleAdjacency(edges);

%Split off border connections
ind = conn(:,2)==0;
border = conn(ind,1);
conn(ind,:) = [];

%Set up sparse connection array and find components <<Graph_Components.m>>
shells = Graph_Components(conn, size(F, 1));

%Only run the following if F is requested as an output
if nargout>=2
    
    %Split F into shells
    F_new = cell(max(shells), 1);
    for ind = 1:numel(F_new)
        F_new{ind} = F(shells==ind,:);
    end
    
    %Modify V if provided, store in output <<DeleteUnreferencedElements.m>>
    if V_in
        [varargout{1}, F_new] = cellfun(@(x) DeleteUnreferencedElements(V, x), ...
                                        F_new, 'UniformOutput', false);
    end
end

%If requested, determine whether or not shells are closed
if nargout-V_in>=3
    border = shells(border);
    closed = true(max(shells), 1);
    closed(unique(border)) = false;
    varargout{nargout-2} = closed;
end


% %Original code (wb20061205), begin:
% 
% %Find triangle connections <<TRI_DetermineTriangleConnections.m>>
% if nargin>=2 && isequal(varargin{2}, 'connections')
%     connections = varargin{1};
% else
%     F = varargin{1};
%     connections = TRI_DetermineTriangleConnections(F);
%     if nargin>=2
%         V = varargin{2};
%     end
% end
% 
% %Initialise while-loop
% shells  = zeros(size(varargin{1}, 1), 1); %Indicates to which shell each triangle belongs
% triangle_start = 1;                       %Indicates triangle belonging to the first shell
% F_new = {};                               %Cell array containing separate F-matrices
% 
% %Loop over all seperate shells
% while ~isempty(triangle_start)
%     
%     %Split off the next shell <<TRI_SplitOffShell.m>>
%     ind = TRI_SplitOffShell(connections, triangle_start);
%     shells(ind) = max(shells)+1;
%     
%     %Find triangles which haven't been put in a shell yet
%     triangle_start = find(shells==0, 1);
%     
%     %Split off part from F and place it in F_new (if F_new is needed)
%     if nargout>1
%         F_new{end+1, 1} = F(ind,:);
%     end
% end
% 
% %Store shells in varargout
% varargout{1} = shells;
% 
% %If only one output argument was requested, return
% if nargout<=1
%     return;
% end
% 
% %If V is provided in the input, modify F_new and create V_new
% if nargin>=2
%     
%     %Keep only vertices belonging to each shell, and re-reference F_new
%     %<<DeleteUnreferencedElements.m>>
%     [V_new, F_new] = cellfun(@(x) DeleteUnreferencedElements(V, x), F_new, 'UniformOutput', false);
%     
%     %Store V_new in varargout
%     varargout{3} = V_new;
% end
% 
% %Store F_new in varargout
% varargout{2} = F_new;
% 
% %If requested, determine whether or not shells are closed
% %<<TRI_DetermineBorderEdges.m>>
% if nargout>=nargin+2
%     varargout{nargout} = cellfun(@(x) isempty(TRI_DetermineBorderEdges(x)), F_new);
% end
% 
% %Original code (wb20061205), end.