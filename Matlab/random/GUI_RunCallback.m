function GUI_RunCallback(src, evt, varargin)
%GUI_RunCallback  Run a callback stored in an object.
%wb20060925
%
%   Syntax:
%    GUI_RunCallback(src, evt, args{:})
%
%   Input:
%    src:  Handle to the object generating the callback. Automatically
%          passed by Matlab when function handle syntax is used.
%    evt:  Event data structur. Automatically passed by Matlab when
%          function handle syntax is used.
%    args: Contents of any callback property.
%
%   Output: none
%
%   Effect: This function will run the provided callback as if it were run
%   by a handle graphics object.
%
%   Dependencies: none
%
%   Known parents: GUI_FreeRotate.m
%                  GUI_VisualisationKeyPress.m
%                  GUI_SelectContours.m

%Created on 25/09/2006 by Ward Bartels.
%Stabile, fully functional.


%Fail-safe
if isempty(varargin) || isempty(varargin{1})
    return
end

%Split up varargin in function and arguments
fun = varargin{1};
args = varargin(2:end);

%Evaluate fun without arguments or create handle from fun string
if ischar(fun)
    if isempty(args)
        evalin('base', fun);
    else
        fun = str2func(fun);
    end
end

%Evaluate function indicated by handle
if isa(fun, 'function_handle')
    fun(src, evt, args{:});
end