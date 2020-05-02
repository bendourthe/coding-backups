function [stiffness, rotax, loc, varargout] = KIN_DirectedSpringRotation(vertices, directions, force_c, sorted)
%KIN_DirectedSpringRotation  Find rotational eigenmodes for spring system.
%
%   Syntax:
%    [stiffness, rotax, loc] = ...
%        KIN_DirectedSpringRotation(vertices, directions, force_c, sorted)
%
%   Input:
%    vertices:   N-by-3 array defining vertices where springs are attached.
%                The rows correspond to different springs and the columns
%                correspond to X-, Y- and Z-coordinates. The elements are
%                coordinate values.
%    directions: N-by-3 array defining normal vectors that indicate the
%                springs' directions. The rows correspond to different
%                springs and the columns correspond to X-, Y- and Z-
%                components. The elements are vector components.
%    force_c:    N-element column vector defining spring force constants.
%                The elements correspond to different springs. Optional,
%                defaults to ones(size(vertices, 1), 1).
%    sorted:     Logical scalar indicating whether or not stiffness must be
%                sorted (in ascending order). Optional, defaults to true.
%
%   Output:
%    stiffness: 3-element column vector containing rotational stiffnesses
%               expressed in torque units per radian. The elements
%               correspond to different eigenmodes.
%    rotax:     3-by-3 array containing unit-length rotation axes. The rows
%               correspond to different eigenmodes and the columns
%               correspond to X-, Y- and Z-components. The elements are
%               vector components expressed in radians per time unit.
%    loc:       3-by-3 array containing points on the rotation screw axes.
%               The rows correspond to different eigenmodes and the columns
%               correspond to X-, Y- and Z-coordinates. The elements are
%               coordinate values.
%
%   Effect: This function will find three rotational eigenmodes in a system
%   consisting of an object bound by a number of directed springs. It will
%   compose a 3-by-3 correlation matrix for the entire system. This matrix
%   will have 3 eigenmodes. The eigenvalues denote the stiffness, and the
%   eigenvectors are used to calculate the rotational screw axis of the
%   corresponding eigenmode.
%
%   References: Lin, Q.; Burdick, J.W.; Rimon, E. A Stiffness-Based Quality
%               Measure for Compliant Grasps and Fixtures. IEEE Trans Rob
%               Autom 16(6), 2000, 675-688.
%
%   Dependencies: KIN_ScrewAxis.m
%
%   Known parents: Femur_ShaftAxis.m


%Created on 31/01/2008 by Ward Bartels.
%WB, 10/03/2009: Removed output of pitch.
%WB, 21/12/2009: Added sorted input argument.
%Stabile, fully functional.


%Calculate moment arms
momarms = cross(vertices, directions, 2);

%Apply weights if necessary
if nargin<3
    momarms_w = momarms;
    directions_w = directions;
else
    weights = force_c(:,[1 1 1]);
    momarms_w = momarms.*weights;
    directions_w = directions.*weights;
end

%Compose correlation matrices
Crr = momarms.'*momarms_w;
Crt = momarms.'*directions_w;
Ctt = directions.'*directions_w;
Ctf = -Ctt\Crt.'; %relates rotation to translation; used later
C = Crr+Crt*Ctf;

%Calculate eigenvalues and eigenvectors
[eigvec, eigval] = eig(C);

%Rearrange eigenvalues and eigenvectors into stiffness and rotation axes
stiffness = diag(eigval, 0);
rotax = eigvec.';

%Sort stiffnesses if requested
if nargin<4 || sorted
    [stiffness, ind] = sort(stiffness);
    rotax = rotax(ind,:);
end

%Calculate location of screw axis <<KIN_ScrewAxis.m>>
if nargout>=3
    loc = KIN_ScrewAxis(rotax, rotax*Ctf.');
end