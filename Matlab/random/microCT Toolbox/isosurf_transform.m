% ########################################################################
% # Name:              isosurf_transform.m (v1.0)                        #
% # Purpose:           Tranforms input isosurface object                 #
% # Author:            Borys Drach                                       #
% # Created:           09/06/12                                          #
% # Copyright:         (c) 2012 Computational Mechanics Lab              #
% #                             Mechanical Engineering Department        #
% #                             University of New Hampshire              #
% ########################################################################

function [ poresurf, sf ] = isosurf_transform( poresurf, V )

% The input 3D shape 'poresurf' object is: a) rotated using rotation matrix 'V' so that principal axes are aligned with the coordinate system axes (I11 is aligned with x-axis, I22 - with y-axis, I33 - with z-axis)
% b) resized so that object's largest dimension is equal to 2 (for FEA analysis)
% c) moved so that object's midpoint lies at the origin of coordinates

for i = 1:size(poresurf.vertices,1)         % Initiation of the loop that goes through all vertices in the 3D object
    poresurf.vertices(i,:) = V'*[poresurf.vertices(i,1);poresurf.vertices(i,2);poresurf.vertices(i,3)];     % rotation of the object
end

% Maximum and minimum coordinates of the object after rotation
max_X = max(max(poresurf.vertices(:,1)));
min_X = min(min(poresurf.vertices(:,1)));
max_Y = max(max(poresurf.vertices(:,2)));
min_Y = min(min(poresurf.vertices(:,2)));
max_Z = max(max(poresurf.vertices(:,3)));
min_Z = min(min(poresurf.vertices(:,3)));

% Half of linear dimensions of the object after rotation
a(1) = (max_X-min_X)/2;
a(2) = (max_Y-min_Y)/2;
a(3) = (max_Z-min_Z)/2;

% Scale factor is determined from the largest linear dimension
sf = max(a);

% All vertices are moved so that the midpoint of the object is located at the origin of coordinates. The object is scaled by the scale factor 'sf' so that largest linear dimension of the object is equal to 2 
poresurf.vertices(:,1) = (poresurf.vertices(:,1) - (min_X+a(1))*ones(size(poresurf.vertices,1),1))/sf;
poresurf.vertices(:,2) = (poresurf.vertices(:,2) - (min_Y+a(2))*ones(size(poresurf.vertices,1),1))/sf;
poresurf.vertices(:,3) = (poresurf.vertices(:,3) - (min_Z+a(3))*ones(size(poresurf.vertices,1),1))/sf;

end