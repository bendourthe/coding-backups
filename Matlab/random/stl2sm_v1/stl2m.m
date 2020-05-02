function filename = stl2m
% Copyright 2011 The MathWorks, Inc.
  
  try
    [filename, pathname] = uigetfile('*.stl');
    %Change filename from stl to m
    
    if exist(fullfile(pathname,[filename(1:end-4) '.m']),'file')
      disp(['using already converted file ' [filename(1:end-4) '.m']]);
      filename = [filename(1:end-4) '.m'];
    else
      fidr = fopen(filename,'r');
      cr = onCleanup(@()fclose(fidr));
      
      T=textscan(fidr,'%s');
      idx=strfind(T{:},'vertex');
      a=zeros(length(idx),1);
      for k = 1:length(idx)
        if ~isempty(idx{k})
          a(k)=1;
        end
      end
      
      idx=find(a);
      filename = [filename(1:end-4) '.m'];
      fidw = fopen(filename,'w+');
      cw = onCleanup(@()fclose(fidw));
      
      fprintf(fidw,'%s\n','k=1;');
      
      T=T{1};
      j=1;
      for k=1:length(idx)/3
        fprintf(fidw,'triangle(k).v1 = [%s %s %s];\n',T{idx(j)+1},T{idx(j)+2},T{idx(j)+3});
        j=j+1;
        fprintf(fidw,'triangle(k).v2 = [%s %s %s];\n',T{idx(j)+1},T{idx(j)+2},T{idx(j)+3});
        j=j+1;
        fprintf(fidw,'triangle(k).v3 = [%s %s %s];\n',T{idx(j)+1},T{idx(j)+2},T{idx(j)+3});
        j=j+1;
        fprintf(fidw,'%s\n','k=k+1;');
      end
    end
  catch me
    me
  end
  
end
