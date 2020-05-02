function angles = VectorAngles(base, vect, rotref)
%VectorAngles  Calculate directed angles between vectors.
%wb20060615
%
%   Syntax:
%    angles = VectorAngles(base, vect, rotref)
%
%   Input:
%    base:   2- or 3-element row vector containing a base vector. The
%            elements represent X-, Y- and possibly Z-coordinates.
%    vect:   N-by-2 or N-by-3 array containing vectors. Each row represents
%            one vector; the columns contain X-, Y- and possibly Z-
%            coordinates.
%    rotref: 3-element row vector containing a rotation reference vector.
%            Will be set to [0 0 1] if size(base, 2)<3. Optional, defaults
%            to cross(base, vect(1,:)).
%
%   Output:
%    angles: Column vector containing the angles in radians between the
%            base vector and the vectors in vect. Element n of angles will
%            be negative if dot(cross(base, vect(n,:)), rotref)<0. The
%            elements of angles lie between pi and -pi.
%
%   Effect: This function will calculate directed angles between a set of
%   vectors and a base vector. The sign of the angles will depend on the
%   direction of the cross product of the base vector and the vectors,
%   relative to the rotation reference.
%
%   Dependencies: none
%
%   Known parents: TRI_RegionCutter.m

%Created on 15/06/2006 by Ward Bartels.
%Stabile, fully functional.


%Add third dimension if necessary
if size(base, 2)<3
    base(:,3) = 0;
    vect(:,3) = 0;
    rotref = [0 0 1];
end

%Calculate rotation reference vector if necessary
if ~exist('rotref', 'var')
    rotref = cross(base, vect(1,:));
end

%Calculate the absolute values of the angles
angles = real(acos(vect*base.'./sqrt(sum(vect.^2, 2))/norm(base)));

%Check which rotation vectors are oriented against the rotation reference
ind = cross(repmat(base, size(vect, 1), 1), vect, 2)*rotref.'<0;
angles(ind) = -angles(ind);