function data=read_data(path,file,ch);
filepath=strcat(path,filesep,file);

A = fopen(filepath, 'r');
B = fread(A,'float');
 
[totalcal, to] = size(B);

m = totalcal/(ch+1);
data = reshape(B, ch+1, m);
data = data';
ST = fclose(A);