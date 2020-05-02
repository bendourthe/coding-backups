% ########################################################################
% # Name:              isosurf_geom_params.m (v1.0)                      #
% # Purpose:           Computes geometrical parameters of input          #
% #                    isosurface object                                 #
% # Author:            Borys Drach                                       #
% # Created:           09/06/12                                          #
% # Copyright:         (c) 2012 Computational Mechanics Lab              #
% #                             Mechanical Engineering Department        #
% #                             University of New Hampshire              #
% ########################################################################

function [dim_Xi,dim_Yi,dim_Zi,atotal,Vol,Xgtotal,J,I11,I22,I33,a,b,c,V] = isosurf_geom_params(poresurf)

% Function calculates:
%     - dimensions of the extracted 3D shape (dim_Xi, dim_Yi, dim_Zi)
%     - object's surface area (atotal), 
%     - volume (Vol), 
%     - center of mass (Xgtotal), 
%     - inertia tensor (J), 
%     - 3 principal moments of inertia (I11, I22, I33),
%     - semiaxes of the approximating ellipsoid (having the same principal moments of inertia and volume as the original object) - a,b,c
%     - axes of the principal moments of inertia (V)

% Number of elements in the object
numelms = size(poresurf.faces,1);

S = 1/120*[2,1,1;1,2,1;1,1,2];                                              % Covariance of the canonical tetrahedron (see paper)

% Memory allocation for element areas, volumes, 
ai = zeros(1,numelms);
Voli = zeros(1,numelms);
Ctot = zeros(3,3);
Xgi = zeros(3,numelms);
Xgtotal = zeros(3,1);

% Dimensions of the extracted 3D shape
dim_Xi = max(poresurf.vertices(:,1)) - min(poresurf.vertices(:,1));
dim_Yi = max(poresurf.vertices(:,2)) - min(poresurf.vertices(:,2));
dim_Zi = max(poresurf.vertices(:,3)) - min(poresurf.vertices(:,3));

% In the loop, every face of the 3D shape is considered. Vectors on which the faces are constructed are determined 
for i=1:numelms
    v1(1,:) = poresurf.vertices(poresurf.faces(i,1),:); 
    v2(1,:) = poresurf.vertices(poresurf.faces(i,2),:);
    v3(1,:) = poresurf.vertices(poresurf.faces(i,3),:);
    
    V(1,:) = v1(1,:); V(2,:) = v3(1,:); V(3,:) = v2(1,:);
    V = V';
    
    ai(i) = 1/2*abs(norm(cross((v2-v1),(v3-v1))));	% Area of the element is determined as half the norm of the cross product of two vectors constructed on the vertices of the element
    Voli(i) = 1/6*det(V);                           % Volume of the element is determined as 1/6 of the mixed product of three vectors constructed on the vertices of the element (see paper)
    Xgi(:,i) = mean([zeros(3,1) V],2);              % Center of mass of the tetrahedron with surface element as base and 4th point at the origin of coordinates
    Ctot = Ctot + det(V)*V*S*V';                    % Covariance of the object (see paper)
    Xgtotal = Xgtotal + Xgi(:,i)*Voli(i);           % Center mass of the object
end

% Total surface area of the object is determined as sum of all element areas comprising the object
atotal = sum(ai);

% Total volume of the object = sum of all element volumes
Vol = sum(Voli);

% Center of mass of the object
Xgtotal = Xgtotal/Vol;

% Covariance of the object calculated around the origin of coordinates (see paper for details)
dx = -Xgtotal;
Ctot = Ctot + Vol*(dx*Xgtotal'+Xgtotal*dx'+dx*dx');

% Inertia tensor calculated based on the covariance tensor (see paper)
J = eye(3)*trace(Ctot) - Ctot;

% Principal moments of inertia and principal axes
A = [J(1,1), J(1,2), J(1,3); J(1,2), J(2,2), J(2,3); J(1,3), J(2,3), J(3,3)];
[V,D] = eig(A);
I11 = D(1,1); I22 = D(2,2); I33 = D(3,3);

% Semiaxes of the approximating ellipsoid (having the same principal moments of inertia and volume as the original object)
a = sqrt(5/2/Vol*(I22+I33-I11));
b = sqrt(5/2/Vol*(I33+I11-I22));
c = sqrt(5/2/Vol*(I11+I22-I33));
end