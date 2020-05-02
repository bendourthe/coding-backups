% ########################################################################
% # Name:              filename_app_zeros.m (v1.0)                       #
% # Purpose:           Zeros in filename                                 #
% # Author:            Borys Drach                                       #
% # Created:           09/06/12                                          #
% # Copyright:         (c) 2012 Computational Mechanics Lab              #
% #                             Mechanical Engineering Department        #
% #                             University of New Hampshire              #
% ########################################################################

function [ zeross ] = filename_app_zeros( num_digs, num )
    % String 'zeross' will be appended to the filename for consistent numbering of object files
    zeross = '';                            % Initialization of the string
    for i = 1:num_digs-size(num2str(num),2) % Loop from 1 to the number of zeros (based on the object number 'num' and total number of digits in the filename requested by the user 'num_digs')
        zeross = strcat(zeross,'0');        % In every cycle of the loop, '0' is added to the 'zeross' string
    end
end