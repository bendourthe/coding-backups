function [totalVolume,totalArea] = stlVolume(V,F)
% Given a surface triangulation, compute the volume enclosed using
% divergence theorem.
% Assumption:Triangle nodes are ordered correctly, i.e.,computed normal is outwards
% Input: V: (N-by-3 array defining vertices), F: (M-by-3 array defining a triangle mesh)
% Output: total volume enclosed, and total area of surface  
% Author: K. Suresh; suresh@engr.wisc.edu

% Compute the vectors d13 and d12
d13= [(V(F(:,2),1)-V(F(:,3),1)), (V(F(:,2),2)-V(F(:,3),2)),  (V(F(:,2),3)-V(F(:,3),3))];
d12= [(V(F(:,1),1)-V(F(:,2),1)), (V(F(:,1),2)-V(F(:,2),2)), (V(F(:,1),3)-V(F(:,2),3))];
cr = cross(d13.',d12.');%cross-product (vectorized)
area = 0.5*sqrt(cr(1,:).^2+cr(2,:).^2+cr(3,:).^2);% Area of each triangle
totalArea = sum(area);
crNorm = sqrt(cr(1,:).^2+cr(2,:).^2+cr(3,:).^2);
zMean = (V(F(:,1),3)+V(F(:,2),3)+V(F(:,3),3))/3;
nz = -cr(3,:)./crNorm;% z component of normal for each triangle
volume = area.*zMean.'.*nz; % contribution of each triangle
totalVolume = sum(volume);%divergence theorem