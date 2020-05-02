function [points] = iv2stack_fast(ivFile,stackFile)
%sorts points in ascending order; 
%edited 2-10-2012

[points coordIndex] = read_vrml_fast(ivFile);

points=sortrows(points,3);
if (exist('stackFile','var')==1),
    dlmwrite(stackFile,points,'delimiter', '\t','newline', 'pc','precision', '%.4f');
end;