function [R, T, Vdat, err] = ICP_RegisterMatchedDatasets(Vdat, Vmod, weights, method)
%ICP_RegisterMatchedDatasets  Register matched point clouds onto eachother.
%
%   Syntax:
%    [R, T, Vdat, err] = ...
%        ICP_RegisterMatchedDatasets(Vdat, Vmod, weights, method)
%
%   Input:
%    Vdat:    M-by-N array defining vertices in the data set. The rows
%             correspond to different data set vertices and the columns
%             correspond to different dimensions. The elements are
%             coordinate values.
%    Vmod:    M-by-N array defining vertices in the model set. The rows
%             correspond to different model set vertices and the columns
%             correspond to different dimensions. The elements are
%             coordinate values.
%    weights: M-element column vector containing weight factors. The
%             elements correspond to different data set and model set
%             vertices. Will be disregarded if empty. Optional, defaults to
%             an empty matrix.
%    method:  String indicating the type of transformation to be used. Can
%             be set to:
%              'rigid':     Rigid transformation.
%              'isometric': Rigid + mirror reflection.
%              'scaling':   Rigid + uniform scaling.
%              'similar':   Rigid + reflection + scaling.
%              'affine':    Rigid + reflection + scaling + shearing.
%             Optional, defaults to 'rigid'.
%
%   Output:
%    R:    N-by-N transformation matrix. The rows correspond to input
%          coordinates and the columns correspond to output coordinates.
%    T:    N-element translation row vector. The elements correspond to
%          coordinates.
%    Vdat: M-by-N array defining vertices in the transformed data set. The
%          rows correspond to different data set vertices and the columns
%          correspond to different dimensions. The elements are coordinate
%          values.
%    err:  Mean squared error (distance between paired points).
%
%   Effect: This function will calculate the transformation matrix R and
%   the translation vector T that register the data set onto the model set.
%   The transformation Vdat*R+repmat(T, size(Vdat, 1), 1) is the
%   least-squares best transformed approximation for Vmod.
%
%   References: Challis, J.H. A procedure for determining rigid body
%               transformation parameters. J Biomech 28(6), 1995, 733-737.
%               Söderkvist, I.; Wedin, P-Å. Determining the movements of
%               the skeleton using well-configured markers. J Biomech
%               26(12), 1993, 1473-1477.
%               Umeyama, S. Least-Squares Estimation of Transformation
%               Parameters Between Two Point Patterns. IEEE Trans Pat Anal
%               Mach Int 13(4), 1991, 376-380.
%
%   Dependencies: SquaredDistanceFromVertexToVertex.m
%
%   Known parents: ICP_RegisterDatasets.m
%                  ICP_RegisterSurfaces.m
%                  KIN_RegisterMarkers.m
%                  Muscle_SelectDatabaseBone.m

%Created on 29/11/2006 by Ward Bartels.
%WB, 03/01/2007: Added weights.
%WB, 16/05/2007: Improved handling of weights.
%WB, 21/05/2007: Added affine and isometric transformations.
%WB, 09/07/2009: Added scaling and similar transformations.
%WB, 13/04/2010: Fixed error in weighting of affine transformations.
%Stabile, fully functional.


%Set default method
if nargin<4, method = 'rigid'; end

%Determine types of transformations required
rigid = strcmpi(method, 'rigid');
isometric = strcmpi(method, 'isometric');
scaling = strcmpi(method, 'scaling');
similar = strcmpi(method, 'similar');
affine = strcmpi(method, 'affine');
noreflect = scaling || rigid;
stretch = scaling || similar;

%Check if weights were provided
if nargin<3 || isempty(weights)
    
    %Define inactive weight application function
    weightfun = @(x) x;
    
    %Define unweighted mean function
    meanfun = @(x) mean(x, 1);
    
else
    
    %Normalise weights
    weights = weights/sum(weights);
    
    %Define weighted mean function
    meanfun = @(x) weights.'*x;
    
    %Modify weights for affine calculation (will not affect meanfun)
    if affine
        weights = sqrt(weights);
    end
    
    %Define weight application function
    weightfun = @(x) x.*weights(:,ones(size(x, 2), 1));
end

%Create column vector with a 1 for each point pair
n = ones(size(Vdat, 1), 1);

%Calculate centroids of point clouds
cendat = meanfun(Vdat);
cenmod = meanfun(Vmod);

%Store cenmod and expand
T = cenmod;
cenmod = n*cenmod;

%Recenter point clouds
Vdat = Vdat-n*cendat;
Vmod = Vmod-cenmod;

%Check if affine transformation was requested
if rigid || isometric || scaling || similar
    
    %Singular value decomposition of weighted covariance matrix
    [U, S, V] = svd(Vmod.'*weightfun(Vdat));
    
    %Compensate reflection if requested
    if noreflect
        U(:,end) = U(:,end)*det(U)*det(V);
    end
    
    %Calculate the rotation matrix
    R = V*U.';
    
    %Apply scale factor if requested
    if stretch
        R = R*(trace(S)/sum(weightfun(sum(Vdat.^2, 2))));
    end
    
elseif affine
    
    %Calculate affine transformation matrix
    R = weightfun(Vdat)\weightfun(Vmod);
    
else
    error([mfilename ':UnknownMethod'], 'Unknown method.');
end

%Calculate translation vector
T = T-cendat*R;

%Transform the dataset if necessary
if nargout>=3
    Vdat = Vdat*R+cenmod;
end

%Calculate the error if necessary <<SquaredDistanceFromVertexToVertex.m>>
if nargout>=4
    err = meanfun(SquaredDistanceFromVertexToVertex(Vdat, Vmod+cenmod, true));
end