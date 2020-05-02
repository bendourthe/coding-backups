function [cop_x,cop_y] = COP_trajectory_Pedar(data, insolearea, subj_cent)
%% DESCRIPTION:
%
%
% [cop_x,cop_y] = COP_trajectory_Pedar(data, insolearea, subj_cent)

% Calculates the COP of a pressure insole (Pedar, Orpyx)
%
%   Author: Anna Meixner
%
%   Last update: Feb. 06th, 2018
%
%% Input:
%   - data: cut data of one pressure insole (left or right) with the sensors in
%           the columns and the time points in the rows
%   - insolearea: area of each sensor as a column vector
%   - subj_cent: 1st column x-coordinate of the center of each sensor,
%                2nd column y-coordinate of the center of each sensor
%   
%% Output:
%   - cop_x: x-coordinate of the center of pressure
%   - cop_y: y-coordinate of the center of pressure
%% force per sensor
    ascholdingarray=zeros(size(data,1),size(data,2));

    for l = 1:size(data,2)
          ascholdingarray(:,l) =data(:,l).*insolearea(l);   
    end

    pedar_fsum = squeeze(sum(ascholdingarray,2));

    %for not deviding through 0
    if find(pedar_fsum ==0)~= 0 
        pedar_fsum((pedar_fsum ==0),1)= 0.000000000000000000000000000000000000000000000000000001; 
    end
    force_matrix= ascholdingarray; %Force matrix for whole insole
%     cent_matx = repmat(subj_cent(1,:),size(ascholdingarray,1),1);
%     cent_maty= repmat(subj_cent(2,:),size(ascholdingarray,1),1);
    cent_matx = repmat(subj_cent(:,1),1,size(ascholdingarray,1))';
    cent_maty = repmat(subj_cent(:,2),1,size(ascholdingarray,1))';
    cop_x = sum((cent_matx.*force_matrix),2)./pedar_fsum; %x-coord COP for whole insole (mm)
    cop_y = sum((cent_maty.*force_matrix),2)./pedar_fsum; %y-coord COP for whole insole (mm)

   %find NAN
   p=0;
   indices = find(isnan(cop_x) == 1);
   if size(indices,1)>0
       for l=1:size(indices,1)
           cop_x (indices(l,1)-p,:)=[];
           cop_y (indices(l,1)-p,:)=[];
           ascholdingarray (indices(l,1)-p,:)=[];
           p=p+1;
       end
   end

      p=0;
   indices = find(isnan(cop_y) == 1);
   if size(indices,1)>0
       for l=1:size(indices,1)
           cop_x (indices(l,1)-p,:)=[];
           cop_y (indices(l,1)-p,:)=[];
           ascholdingarray(indices(l,1)-p,:)=[];
           p=p+1;
       end
   end
end

