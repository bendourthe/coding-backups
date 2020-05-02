function extents = GUI_CoordinateExtents(h)
%GUI_CoordinateExtents  Determine coordinate extents of a group of objects.
%wb20070911
%
%   Syntax:
%    extents = GUI_CoordinateExtents(h)
%
%   Input:
%    h: Column vector of object handles.
%
%   Output:
%    extents: 3-by-2 array containing the coordinate extents of the plot
%             objects under h. The first row contains minima, the second
%             row contains maxima and the three columns correspond to X-,
%             Y- and Z-coordinates. Will contain NaN if no values were
%             found.
%
%   Effect: This function will find all objects under h (including h
%   themselves) that have coordinate data, and return the minimum and
%   maximum coordinates among these.
%
%   Dependencies: none
%
%   Known parents: GUI_ManipulateObjects.m

%Created on 11/09/2007 by Ward Bartels.
%Stabile, fully functional.


%Get handles to all objects under h
h = findobj(h);

%Restrict type of object
h = h(ismember(get(h, 'type'), {'line'; 'patch'; 'surface'}));

%Obtain coordinate values
coords = get(h, {'XData' 'YData' 'ZData'});
coords = cellfun(@(x) x(:), coords, 'UniformOutput', false);
x = vertcat(NaN, coords{:,1});
y = vertcat(NaN, coords{:,2});
z = vertcat(NaN, coords{:,3});

%Get extrema
extents = [min(x) min(y) min(z); max(x) max(y) max(z)];