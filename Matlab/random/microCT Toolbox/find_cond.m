% ########################################################################
% # Name:              find_cond.m (v1.0)                                #
% # Purpose:           Converts 3D matrix element numbers to coordinates #
% # Author:            Borys Drach                                       #
% # Created:           09/06/12                                          #
% # Copyright:         (c) 2012 Computational Mechanics Lab              #
% #                             Mechanical Engineering Department        #
% #                             University of New Hampshire              #
% ########################################################################

function [xx,yy,zz,brdr] = find_cond(V,dims)

max1 = dims(1);
max2 = dims(2);
max3 = dims(3);

szV = size(V,1);

xx = zeros(szV,1);
yy = zeros(szV,1);
zz = zeros(szV,1);

brdr = 0;

for n = 1:szV 
    
    k = floor(V(n)/(max1*max2))+1;
    if (V(n) - (k-1)*max1*max2) == 0
        k = k-1;
        j = max2;
        i = max1;
    else
        j = floor((V(n) - (k-1)*max1*max2)/max1)+1;
        if (V(n)-(k-1)*max1*max2-(j-1)*max1) == 0
           j = j-1;
           i = max1;
        else
           i = V(n) - (k-1)*max1*max2 - (j-1)*max1;
        end
    end
        
    if i == 1 || i == max1 || j == 1 || j == max2 || k == 1 || k == max3
       brdr = 1;
    end
    
    xx(n,1) = j;
    yy(n,1) = i;
    zz(n,1) = k;
    
end

end

