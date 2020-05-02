function d = density()
% Copyright 2011 The MathWorks, Inc.

% d is a map container
% use d('iron') to get iron's density
% density in kg/m^3
  
  material = {'iron','plastic'};
  density  = {  7874, 1000};
  
  d = containers.Map(material, density);
