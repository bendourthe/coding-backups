function array = RunLengthDecode(lengths, values)
%RunLengthDecode  Run-length decoding.
%wb20070723
%
%   Syntax:
%    array = RunLengthDecode(lengths, values)
%
%   Input:
%    lengths: Column vector indicating how often each corresponding element
%             in values needs to be repeated. This variable may contain
%             zeros.
%    values:  N-dimensional array containing the elements used to build the
%             output array. Optional, defaults to (1:numel(lengths)).'.
%
%   Output:
%    array: Column vector obtained by repeating each element from values as
%           often as called for in lengths.
%
%   Effect: This function will perform run-length decoding: each element
%   in lengths determines how many times the corresponding element in
%   values will be repeated to build array. If values is multi-dimensional,
%   its "rows" (e.g. values(i,:,:) if values is 3-D) will be concatenated
%   and replicated in the first dimension. All array types are supported
%   for values.
%
%   Dependencies: none
%
%   Known parents: Contour_SplitInner.m
%                  StackEqualElementIndices.m
%                  Graph_Components.m
%                  BoneModel/BoneModel.m
%                  BoneModel/private/Transform.m
%                  SIMM_ReadBone.m
%                  SIMM_ParseKeywords.m
%                  SIMM_ReadMuscle.m
%                  JointModel/JointModel.m
%                  TRI_CutWithMultiPlane.m
%                  PRI_SphereSegment.m

%Created on 18/09/2006 by Ward Bartels.
%WB, 19/09/2006: Added fail-safe to deal with empty lengths.
%WB, 23/07/2007: Added handling of multi-dimensional values array.
%Stabile, fully functional.


%Remove trailing zeros from lengths
lengths = lengths(1:find(lengths>0, 1, 'last'));

%Fail-safe
if isempty(lengths)
    array = zeros(0, 1);
else
    
    %Calculate cumulative sum of lengths
    lengths_sum = cumsum(lengths);
    
    %Create array with ones where values change
    array = accumarray([1; lengths_sum(1:end-1)+1], ones(size(lengths)), [lengths_sum(end) 1]);
    
    %Build output array
    array = cumsum(array);
end

%Index into values if necessary
if nargin>=2
    subs = repmat({':'}, ndims(values)-1, 1);
    array = values(array, subs{:});
end