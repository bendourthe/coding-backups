function [R, T, Vdat, err, err_rel, iter] = ICP_RegisterDatasets(Vdat, Vmod, Wdat, Wmod, wbar, tes, conn, method, tol, imax)
%ICP_RegisterDatasets  Register point clouds using the ICP algorithm.
%
%   Syntax:
%    [R, T, Vdat, err, err_rel, iter] = ...
%        ICP_RegisterDatasets(Vdat, Vmod, Wdat, Wmod, wbar, tes, conn, ...
%                             method, tol, imax)
%
%   Input:
%    Vdat:   M-by-N array containing vertex coordinates. Each row
%            represents a vertex in the data set; the columns contain
%            coordinate data in different dimensions.
%    Vmod:   M-by-N array containing vertex coordinates. Each row
%            represents a vertex in the model set; the columns contain
%            coordinate data in different dimensions.
%    Wdat:   Column vector containing weight factors for each data point.
%            Will be disregarded if empty. Optional, defaults to an empty
%            matrix.
%    Wmod:   Column vector containing weight factors for each model point.
%            Will be disregarded if empty. Optional, defaults to an empty
%            matrix.
%    wbar:   Logical indicating whether or not a wait bar should be shown.
%            Optional, defaults to true.
%    tes:    Logical indicating whether or not tessellation should be used
%            to speed up the point matching step. Optional, defaults to
%            true.
%    conn:   Scalar KDTreeSearcher object defining a kd-tree over the
%            vertices in Vmod. This variable may be generated using
%            ICP_MatchDatasets(Vmod). It is disregarded if tes is set to
%            false. Optional, will be calculated automatically if empty or
%            omitted.
%    method: String indicating the type of transformation to be used. Will
%            be used as input for ICP_RegisterMatchedDatasets.m; see that
%            function's help text for possible transformations. Optional,
%            defaults to 'rigid'.
%    tol:    Relative tolerance for the mean squared error. If the error
%            changes by a value smaller than tol after one iteration, the
%            registration procedure will stop. Optional, defaults to 1e-4.
%    imax:   Maximum number of iterations. Optional, defaults to 100.
%
%   Output:
%    R:       Square rotation matrix.
%    T:       Translation row vector.
%    Vdat:    M-by-N array containing vertex coordinates. Each row
%             represents a vertex in the registered data set; the columns
%             contain coordinate data in different dimensions.
%    err:     Mean squared error (squared distance between paired points).
%    err_rel: Relative error (change in error through last iteration). A
%             negative value indicates instability.
%    iter:    Number of iterations needed.
%
%   Effect: This function will run the ICP-algorithm in order to find the
%   rotation matrix R and the translation matrix T that register the data
%   set onto the model set. The transformed dataset Vdat*R+T should lie
%   very closely to the modelset Vmod.
%
%   References: Besl, P.J.; McKay, N.D. A Method for Registration of 3-D
%               Shapes. IEEE Trans Pat Anal and Mach Intel 14(2), 1992,
%               239-256.
%
%   Dependencies: ICP_MatchDatasets.m
%                 ICP_RegisterMatchedDatasets.m
%
%   Known parents: Muscle_AlignRegion.m
%                  Muscle_SelectDatabaseBone.m

%Created on 29/11/2006 by Ward Bartels.
%WB, 01/12/2006: Added output of Vdat; waitbar is now immediately shown;
%                relative rather than absolute error is now used.
%WB, 15/12/2006: Speed-up by using tsrchnmx instead of ICP_MatchDatasets.m.
%WB, 03/01/2007: Added weights for data set.
%WB, 01/03/2007: Turned off warning on log of zero.
%WB, 16/05/2007: Added weights for model set.
%WB, 22/05/2007: Added affine and isometric transformations.
%WB, 23/05/2007: Added infinite loop fail-safe and output of err_rel.
%WB, 03/09/2007: Added input of conn.
%WB, 13/04/2010: Same weights are now used in both terms of relative error.
%WB, 05/12/2011: Added call to DRAWNOW after wait bar is closed.
%Stabile, fully functional.


%Set defaults
if nargin<3, Wdat = []; end
if nargin<4, Wmod = []; end
if nargin<5, wbar = true; end
if nargin<6, tes = true; end
if nargin<7, conn = {}; end
if nargin<8, method = 'rigid'; end
if nargin<9, tol = 1e-4; end
if nargin<10, imax = 100; end

%Fail-safe, prevent infinite loop
if tol<0 && imax==Inf
    error([mfilename ':InfiniteLoop'], 'Negative tolerance and infinite iteration cap would cause infinite loop');
end

%Initialise wait bar if requested
if wbar
    wstate = warning('off', 'MATLAB:log:logOfZero');
    wb = waitbar(0, 'Initialising registration...', 'Name', 'Please wait...');
end

%Define anonymous weighted mean and weight application functions
if isempty(Wdat) && isempty(Wmod)
    weightfun = @(pm) [];
    meanfun = @(x, w) mean(x, 1);
else
    if isempty(Wdat)
        weightfun = @(pm) Wmod(pm,1);
    elseif isempty(Wmod)
        weightfun = @(pm) Wdat;
    else
        weightfun = @(pm) Wdat.*Wmod(pm,1);
    end
    meanfun = @(x, w) w.'*x/sum(w);
end

%Store initial data set
Vdat_init = Vdat;

%Find closest point pairs for initial data set
if tes
    
    %Create tessellation and look for closest points
    if isempty(conn)
        conn = ICP_MatchDatasets(Vmod);
    end
    [points_mod, err] = ICP_MatchDatasets(Vdat, Vmod, conn);
    
else
    
    %Exhaustively search for closest points <<ICP_MatchDatasets.m>>
    [points_mod, err] = ICP_MatchDatasets(Vdat, Vmod);
end

%Calculate initial weights and relative error
weights = weightfun(points_mod);
err = meanfun(err, weights);
err_rel = err;

%Initialise iteration counter and transformation
iter = 0;
R = eye(size(Vmod, 2));
T = zeros(1, size(Vmod, 2));

%Update wait bar text
if wbar
    err_i = log(err);
    err_d = err_i-log(tol);
    waitbar(0, wb, 'Running registration - Iteration 1. Error evolution:');
end

%Loop until tolerance or maximum number of iterations is reached
while err_rel>tol && iter<imax
    
    %Refresh graphics, useful if run inside larger loop
    drawnow;
    
    %Increment iteration counter
    iter = iter+1;
    
    %Calculate new transformation <<ICP_RegisterMatchedDatasets.m>>
    [R, T, Vdat] = ICP_RegisterMatchedDatasets(Vdat_init, Vmod(points_mod,:), weights, method);
    
    %Store previous absolute error
    err_rel = err;
    
    %Re-pair points and calculate error <<ICP_MatchDatasets.m>>
    if tes
        [points_mod, err] = ICP_MatchDatasets(Vdat, Vmod, conn, points_mod);
    else
        [points_mod, err] = ICP_MatchDatasets(Vdat, Vmod);
    end
    
    %Calculate relative error using previous weights
    err_rel = err_rel-meanfun(err, weights);
    
    %Calculate new weights and absolute error
    weights = weightfun(points_mod);
    err = meanfun(err, weights);
    
    %Update wait bar if necessary
    if wbar
        waitbar((err_i-log(err_rel))/err_d, wb, ...
                ['Running registration - Iteration ' num2str(iter+1) '. Error evolution:']);
    end
end

%Close wait bar if necessary
if wbar
    warning(wstate);
    close(wb);
    drawnow('expose');
end