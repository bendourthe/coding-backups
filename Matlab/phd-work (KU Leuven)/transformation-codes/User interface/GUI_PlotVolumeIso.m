function [obj, ax, thresh] = GUI_PlotVolumeIso(x, y, z, vol, numsurf, colors)
%GUI_PlotVolumeIso  Plot a volumetric dataset with isosurfaces.
%
%   Syntax:
%    [obj, ax, thresh] = GUI_PlotVolumeIso(x, y, z, vol, numsurf, color)
%    [obj, ax, thresh] = GUI_PlotVolumeIso(x, y, z, vol, numsurf, colors)
%
%   Input:
%    x:       L x 1 double array defining coordinate values on X-axis.
%              Rows:     X-coordinates
%              Contents: Coordinate values
%    y:       M x 1 double array defining coordinate values on Y-axis.
%              Rows:     Y-coordinates
%              Contents: Coordinate values
%    z:       N x 1 double array defining coordinate values on Z-axis.
%              Rows:     Z-coordinates
%              Contents: Coordinate values
%    vol:     L x M x N double array defining volumetric scalar data.
%              1st dim:  X-coordinates
%              2nd dim:  Y-coordinates
%              3rd dim:  Z-coordinates
%              Contents: Values of scalar field over X, Y and Z
%    numsurf: Scalar double defining the number of isosurfaces to be shown.
%    color:   1 x 3 double array defining color for all isosurfaces.
%              Columns:  Red, green and blue
%              Contents: Color values, between 0 and 1
%    colors:  numsurf x 3 double array defining different colors for the
%             isosurfaces.
%              Rows:     Isosurfaces that will be shown
%              Columns:  Red, green and blue
%              Contents: Color values, between 0 and 1
%
%   Output:
%    obj:    numsurf x 2 double array containing patch handles.
%             Rows:     Isosurfaces that are shown
%             Columns:  Surfaces and bordering lines
%             Contents: Handles to patch objects displaying isosurfaces and
%                       bordering lines
%    ax:     Scalar double containing a handle to the axes on which the
%            plot is shown.
%    thresh: numsurf x 1 double array containing threshold values used for
%            isosurfaces.
%             Rows:     Isosurfaces that are shown
%             Contents: Threshold values used to calculate isosurfaces, in
%                       the same units as vol.
%
%   Effect: This function will visualize a volumetric dataset, i.e. a
%   scalar field in 3 dimensions defined by voxel values. A number of
%   isosurfaces, specified in numsurf, is calculated. Along these surfaces,
%   vol has the same value. These isosurfaces are visualized transparently
%   with the specified color. Also, the borders of the isosurfaces are
%   calculated and shown as black lines.
%
%   Dependencies: TRI_BorderEdges.m
%                 GUI_PlotShells.m
%
%   Known parents: none

%Created on 05/11/2007 by Ward Bartels.
%Stabile, fully functional.


%Get threshold values
thresh = linspace(min(vol(:)), max(vol(:)), numsurf+2).';
thresh([1 end]) = [];

%Calculate all isosurfaces
[F, V] = arrayfun(@(t) isosurface(y, x, z, vol, t), thresh, 'UniformOutput', false);
V = cellfun(@(v) v(:,[2 1 3]), V, 'UniformOutput', false);

%Get isosurface edges <<TRI_BorderEdges.m>>
edges = cellfun(@TRI_DetermineBorderEdges, F, 'UniformOutput', false);

%Plot isosurfaces with edges <<GUI_PlotShells.m>>
obj = zeros(numsurf, 2);
[obj(:,1), ignoble, ax] = GUI_PlotShells(F, V, colors);
set(obj(:,1), 'FaceAlpha', 0.55);
obj(:,2) = cellfun(@(e, v) patch('Vertices', v, 'Faces', e, 'Parent', ax, 'FaceAlpha', 0, 'EdgeAlpha', 1), edges, V);

% %Adjust normals
% normals = arrayfun(@(o) isonormals(x, y, z, vol, o), obj, 'UniformOutput', false);
% set(obj, {'VertexNormals'}, normals);