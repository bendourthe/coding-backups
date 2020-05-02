function answer = GUI_QuestionDialog(question, buttons, title, default, modal)
%GUI_QuestionDialog  Modal or non-modal question dialog box.
%wb20060324
%
%   Syntax:
%    answer = GUI_QuestionDialog(question, buttons, title, default, modal)
%
%   Input:
%    question: String containing the question which will appear in the
%              dialog box.
%    buttons:  N-by-1 cell array of strings containing the button names. N
%              may range from 1 to 3. Optional, defaults to {'Yes'; ...
%              'No'; 'Cancel'}.
%    title:    String containing the title of the dialog box. Optional,
%              defaults to an empty string.
%    default:  Integer indicating the location of the default selection in
%              buttons. This variable may range from 1 to length(buttons).
%              Optional, defaults to 1.
%    modal:    Logical indicating whether or not the dialog box should be
%              modal (i.e. the user is prevented from manipulating any
%              other window while the dialog box is open). Optional,
%              defaults to false.
%
%   Output:
%    answer: Integer indicating the location of the selected button in
%            buttons. If the user closed the dialog box, answer will be 0.
%
%   Effect: This function acts as a wrapper for questdlg and GUI_questdlg.
%   It will pop up a dialog box prompting the user to click on one of the
%   buttons. An index to the selected button is returned.
%
%   Dependencies: GUI_questdlg.m
%
%   Known parents: none

%Created on 24/03/2006 by Ward Bartels.
%Stabile, fully functional.


%Set input defaults
if nargin<2, buttons = {'Yes'; 'No'; 'Cancel'}; end
if nargin<3, title = ''; end
if nargin<4, default = 1; end
if nargin<5, modal = false; end

%Pop up dialog box <<GUI_questdlg.m>>
if modal
    name = questdlg(question, title, buttons{:}, buttons{default});
else
    name = GUI_questdlg(question, title, buttons{:}, buttons{default});
end

%Get index into buttons
answer = strmatch(name, buttons, 'exact');
if isempty(answer)
    answer = 0;
else
    answer = answer(1);
end