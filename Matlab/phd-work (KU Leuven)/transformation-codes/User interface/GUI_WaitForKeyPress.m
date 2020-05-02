function GUI_WaitForKeyPress(fig, closefig, keepopen)
%GUI_WaitForKeyPress  Wait for a keypress in a group of figures.
%wb20060313
%
%   Syntax:
%    GUI_WaitForKeyPress(fig, closefig, keepopen)
%
%   Input:
%    fig:      Column vector of figure handles.
%    closefig: Logical, determines whether or not the figures must be
%              closed when a key is pressed in one of them. Optional,
%              defaults to false.
%    keepopen: Logical, determines whether or not the figures must be
%              protected from manual closing by the user as long as no key
%              has been pressed in one of them. Optional, defaults to
%              false.
%
%   Output: none
%
%   Effect: This function will stop the execution of the calling function
%   until the user presses any key in one of the figures indicated by
%   fig. Optionally, the figure windows can be protected from manual
%   closing by the user. Also, it is possible to automatically close all
%   the figure windows as soon as a key is pressed. This function will
%   replace the UserData and CloseRequestFcn properties of the figures.
%
%   Dependencies: none
%
%   Known parents: none

%Created on 13/03/2006 by Ward Bartels.
%Stabile, fully functional.


%Parse input
if nargin==0, fig = gcf; end
closefig = nargin>=2 && closefig;
keepopen = nargin>=3 && keepopen;

%Put handles in UserData of all figures
set(fig, 'UserData', fig);

%Set key press callbak for all figures
set(fig, 'KeyPressFcn', {@KeyPress, closefig});

%Make sure figures stay open if necessary
if keepopen
    set(fig, 'CloseRequestFcn', '');
else
    set(fig, 'CloseRequestFcn', {@KeyPress, true});
end

%Suspend execution until first figure is closed or uiresume is called
uiwait(fig(1));



%Callback, executed when a key is pressed in a figure
function KeyPress(src, evt, closefig)

%Get handles to all figures
fig = get(src, 'UserData');

%Close figures if necessary, otherwise just resume execution
if closefig
    delete(fig);
else
    set(fig, 'CloseRequestFcn', 'closereq');
    set(fig(1), 'waitstatus', 'inactive'); %Alternative to uiresume
end