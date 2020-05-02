function [F_new, V_new, shells] = TRI_CutWithContour(F, V, contour)
%TRI_CutWithContour  Cut mesh along a previously defined contour.
%wb20080104
%
%   Syntax:
%    [F_new, V_new] = TRI_CutWithContour(F, V, contour)
%
%   Input:
%    F:       N-by-3 array containing indices into V. Each row represents a
%             triangle, each element is a link to a vertex in V.
%    V:       N-by-3 array containing vertex coordinates. Each row
%             represents a vertex; the first, second and third columns
%             represent X-, Y- and Z-coordinates respectively.
%    contour: N-by-2 array containing the intersecting contour. Each row
%             represents an edge, each element is a row index into V.
%
%   Output:
%    F_new:  N-by-1 cell array containing separate F-matrices for each
%            shell.
%    V_new:  N-by-1 cell array containing separate V-matrices for each
%            shell.
%
%   Effect: This function will separate a mesh into different shells by
%   cutting through the edges defined by contour. The contour is made up of
%   edges already existing in the mesh.
%
%   Dependencies: TRI_RemoveInvalidTriangles.m
%                 TRI_RemoveBadlyConnectedTriangles.m
%                 TRI_Edges.m
%                 TRI_SeparateShells.m
%
%   Known parents: TRI_RegionCutter.m

%Created on 20/06/2006 by Ward Bartels.
%WB, 05/12/2006: Rewrite based on MatlabBGL.
%WB, 13/02/2007: Added check for badly connected triangles.
%WB, 04/01/2008: Moved code to TRI_SeparateShells.m.
%Stabile, fully functional.


%Remove redundant vertices
[ignoble, ind, jnd] = unique(single(V), 'rows'); %Single precision for thresholding
V = V(ind,:);
jnd = jnd.';
contour = jnd(contour); F = jnd(F);

%Remove invalid triangles <<TRI_RemoveInvalidTriangles.m>>
F = TRI_RemoveInvalidTriangles(F);
[ignoble, ind] = unique(sort(F, 2), 'rows');
F = F(ind,:);

%Remove badly connected triangles <<TRI_RemoveBadlyConnectedTriangles.m>>
triangles = TRI_RemoveBadlyConnectedTriangles(F, V);
F(triangles,:) = [];

%Determine edges <<TRI_Edges.m>>
edges = TRI_Edges(F);

%Eliminate edges on the contour
edges = sort(edges, 2);
edges(ismember(edges, sort(contour, 2), 'rows'),:) = NaN;

%Split up the mesh <<TRI_SeparateShells.m>>
[shells, F_new, V_new] = TRI_SeparateShells(F, V, edges);


% %Original code (wb20061205), begin:
% 
% %Remove redundant vertices
% [ignoble, ind, jnd] = unique(single(V), 'rows'); %Single precision for thresholding
% V = V(ind,:);
% jnd = jnd.';
% contour = jnd(contour); F = jnd(F);
% 
% %Remove invalid triangles <<TRI_RemoveInvalidTriangles.m>>
% F = TRI_RemoveInvalidTriangles(F);
% 
% %Determine edges <<TRI_Edges.m>>
% [edges, triangles] = TRI_Edges(F);
% 
% %Eliminate edges on the contour
% edges = sort(edges, 2);
% edges(ismember(edges, sort(contour, 2), 'rows'),:) = NaN;
% 
% %Add filler element to triangles
% triangles(end+1) = 0;
% 
% %See which edges are common (use reference to filler)
% %<<StackEqualElementIndices.m>>
% stack = StackEqualElementIndices(edges, length(triangles));
% if size(stack, 2)==1
%     stack = [stack ones(size(stack))*length(triangles)]; %Handle mesh containing only border triangles
% end
% connections = triangles(stack(:,[1 2]));
% 
% %Expand connections
% connections = [connections(connections(:,2)~=0,:); fliplr(connections)];
% 
% %Handle faulty meshes
% faulty = nonzeros(triangles(stack(:,3:end)));
% connections = [connections; zeros(size(faulty)) faulty];
% 
% %Sort connections
% [ignoble, ind_sort] = sort(connections(:,2));
% connections = connections(ind_sort,1);
% 
% %Rearrange connections so each row corresponds with a row in F
% connections = reshape(connections, 3, [])';
% 
% %Get shell indices from modified connections <<TRI_SeparateShells.m>>
% shells = TRI_SeparateShells(connections, 'connections');
% 
% %Separate F into shells
% for ind = 1:max(shells)
%     F_new{ind,1} = F(shells==ind,:);
% end
% 
% %Keep only vertices belonging to each shell, and re-reference F_new
% %<<DeleteUnreferencedElements.m>>
% [V_new, F_new] = cellfun(@(x) DeleteUnreferencedElements(V, x), F_new, 'UniformOutput', false);
% 
% % %Original code (wb20061205), end.