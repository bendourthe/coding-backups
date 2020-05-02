function [R, rotax, norm1, norm2] = VectorAligningRotation(vect1, vect2)
%VectorAligningRotation  Create rotation matrix that aligns two vectors.
%
%   Syntax:
%    [R, rotax, norm1, norm2] = VectorAligningRotation(vect1, vect2)
%
%   Input:
%    vect1: 3-element row vector containing the X-, Y- and Z-coordinates of
%           a vector.
%    vect2: 3-element row vector containing the X-, Y- and Z-coordinates of
%           a vector.
%
%   Output:
%    R:     3-by-3 orthonormal rotation matrix.
%    rotax: 3-element row vector containing the X-, Y- and Z-coordinates of
%           the rotation axis of R, which is perpendicular to both vect1
%           and vect2.
%    norm1: 3-element row vector containing the X-, Y- and Z-coordinates of
%           a vector perpendicular to both rotax and vect1.
%    norm2: 3-element row vector containing the X-, Y- and Z-coordinates of
%           a vector perpendicular to both rotax and vect2.
%
%   Effect: This function will create a 3-by-3 orthonormal rotation matrix
%   that transforms vect1 onto vect2, i.e. vect1*R = vect2.
%
%   Dependencies: none
%
%   Known parents: Muscle_AlignRegion.m
%                  PRI_SphereSegment.m

%Created on 04/09/2007 by Ward Bartels.
%Stabile, fully functional.


%Normalize both vectors
vect1 = vect1/norm(vect1);
vect2 = vect2/norm(vect2);

%Find rotation axis that rotates vect1 onto vect2
rotax = cross(vect1, vect2);
if all(abs(rotax)<eps) %vect1 and vect2 are (nearly) parallel
    [ignoble, ind] = min(abs(vect1));
    rotax = cross(vect1, [1 2 3]==ind);
end
rotax = rotax/norm(rotax);

%Create orthonormal rotation matrix
norm1 = cross(vect1, rotax);
norm2 = cross(vect2, rotax);
R = [vect1; rotax; norm1].'*[vect2; rotax; norm2]; %transpose = inverse for unitary matrix