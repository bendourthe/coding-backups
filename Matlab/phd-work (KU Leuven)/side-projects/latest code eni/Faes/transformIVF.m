function [] = transformIV(ivFile, RT, outFile)
% function [] = transformIV(ivFile, RT, outFile)

if (exist('outFile','var')~=1),
    outFile = ivFile;
end;

[pts conn] = read_vrml_fast(ivFile);
conn(:,4) = [];
conn = conn+1;

pts = transformShell(pts,RT,1,1);
% a=pts(:,1);
% b=pts(:,2);
% pts=[-a -b pts(:,3)];
patch2iv(pts,conn,outFile);