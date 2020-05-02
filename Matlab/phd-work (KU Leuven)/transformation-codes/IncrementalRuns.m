function vector = IncrementalRuns(lengths)
%IncrementalRuns  Build vector of incremental runs with specified lengths.
%wb20060920
%
%   Input:
%    lengths: Column vector indicating how long each incremental run should
%             be.
%
%   Output:
%    vector: Column vector built by vertically concatenating
%            (1:lengths(ind)).' for each ind.
%
%   Effect: This function will build up a column vector consisting of
%   multiple incremental runs. Each run is a series of numbers starting at
%   1 and increasing by 1 until the number specified in lengths is reached.
%
%   Dependencies: none
%
%   Known parents: StackEqualElementIndices.m
%                  SIMM_ReadMuscle.m

%Created on 20/09/2006 by Ward Bartels.
%Stabile, fully functional.


%Eliminate zeros from lengths
lengths(lengths==0) = [];

%Fail-safe
if isempty(lengths)
    vector = lengths;
    return
end

%Calculate cumulative sum of lengths
lengths_sum = cumsum(lengths);

%Create array filled with ones, and with jumps where runs start
vector = ones(lengths_sum(end), 1);
vector(lengths_sum(1:end-1)+1) = 1-lengths(1:end-1);

%Build vector
vector = cumsum(vector);