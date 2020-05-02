function varargout = GUI_LoadPNGButtons(tbcolor, varargin)
%GUI_LoadPNGButtons  Load toolbar button icon CData from PNG image.
%wb20060422
%
%   Syntax:
%    [CData1, CData2, ...] = GUI_LoadPNGButtons(tbcolor, file1, file2, ...)
%
%   Input:
%    tbcolor: 3-element row vector containing a color definition describing
%             the background color of the toolbar. The first, second and
%             third columns represent red, green and blue values
%             respectively. Each element is an integer between 0 and 255.
%             If set to an empty matrix, tbcolor will assume the default
%             value of [224 223 227].
%    file#:   String containing a path to a PNG image to be used as toolbar
%             icon.
%
%   Output:
%    CData#: Array containing a toolbar icon in a form that can be used for
%            the "CData" property of uipushtool and uitoggletool objects.
%
%   Effect: This function will read toolbar icons from the set of PNG-files
%   passed as arguments. If the files contain transparency, the toolbar
%   background color defined by tbcolor is automatically blended in. Pixels
%   with the same color as tbcolor are made transparent.
%
%   Dependencies: none
%
%   Known parents: GUI_SimulateSurgery.m
%                  GUI_VisualisationUI.m
%
%   Example:
%    uipushtool('CData', GUI_LoadPNGButtons([224 223 227], 'Cut.png'))

%Created on 22/04/2006 by Ward Bartels.
%Stabile, fully functional.


%Set tbcolor to default if it is empty
if isempty(tbcolor)
    tbcolor = [224 223 227];
end

%Rescale tbcolor to 0-1
tbcolor = tbcolor/255;

%Loop over all button images
varargout = cell(size(varargin));
for ind = 1:nargin-1
    
    %Read image, blend in toolbar background and convert to double
    [image, map] = imread(varargin{ind}, 'PNG', 'BackgroundColor', tbcolor);
    
    %Apply map if necessary, otherwise rescale to 0-1
    if ~isempty(map)
        N = size(map, 1);
        map(1,:) = NaN;
        image = uint16(image)+1;
        image(image==0) = N;
        image = map(cat(3, image, image+N, image+2*N));
    else
        image = double(image)/255;
    end
    
    %Set all transparent pixels to NaN
    image(repmat(all(image==tbcolor(cumsum(ones(size(image)), 3)), 3), [1 1 3])) = NaN;
    
    %Store image in output
    varargout{ind} = image;
end