% ########################################################################
% # Name:              process_vtk_stl.m (v1.0)                          #
% # Purpose:           Processing of 8bit VTK files                    #
% # Author:            Borys Drach                                       #
% # Created:           09/06/12                                          #
% # Copyright:         (c) 2012 Computational Mechanics Lab              #
% #                             Mechanical Engineering Department        #
% #                             University of New Hampshire              #
% ########################################################################

% This code processes black and white VTK file obtained by Computed
% Tomography, finds individual objects, computes their geometrical 
% characteristics and exports meshes of the shapes as STL files

clc;                                    % clears the command window
clear;                                  % clears the workspace
close all;                              % closes all open figure windows

% --------------- INPUT ----------------

init_filename = 'sample_pores';         % name of the VTK file to be processed
ext_stack = '.vtk';                     % VTK file extension
ext_stl = '.stl';                       % extension of the output file. STL - stereolitography format, stores a list of 3-node triangular elements that comprise an object's surface mesh
ext_mat = '.mat';                       % extension of the outpit file that stores geometry information of individual objects

vtk_filename = strcat('vtk\',init_filename, ext_stack);

% Check if VTK file with user specified filename exists
if exist(vtk_filename,'file') == 0
    disp('No such file');
    return
end

II = readvtk(vtk_filename)/255;         % writes (and converts to binary) the information from the VTK file into 3D matrix

dims = size(II);                        % dimensions of the input 3D matrix

cc = bwconncomp(II,6);                  % identifies number of individual objects and stores their voxel coordinates in the 3D matrix (muCT data)
%save('full_setup.mat','-v7.3');

first_num = 168;                        % start with object #
last_num = first_num;                   % finish with object #. 'cc.NumObjects' can be entered here - the total number of objects in the VTK file (see below)
%last_num = cc.NumObjects;              % 
ttl_num = last_num - first_num + 1;     % total number of objects to be analyzed

pix_min = 1;                            % minimum number of elements in the object. If less than this value - ignore the object and move to the next one
pix_max = 1000;                         % maximum number of elements in the object. If more than this value - ignore the object and move to the next one

% default plot settings definition
set(0,'DefaultFigureColor','white', 'DefaultFigurePaperUnits', 'inches',...
    'DefaultFigureResize', 'off','DefaultFigurePaperPositionMode', 'manual',...
    'DefaultAxesPosition',[0 0 1 1]);

ii = 0;         % counter, number of successfully written objects
oo = 0;         % counter, number of objects with dimension less than 1px in any of the three axes
pp = 0;         % counter, number of objects with dimension less than 1px in any of the three axes
m = 1;          % integer loop variables

dat_table = zeros(ttl_num,26);         % memory allocation for the table where pore paramteres will be stored

% --------------- PROCESSING ---------------- 
for num = first_num:last_num              % start of the loop through the range of objects to be analyzed.

    % In this program objects that touch surfaces of the imported 3D matrix (contatining all the objects to be analyzed) are eleminated from the analysis
    % Since 'bwconncomp' function outputs element numbers of the connected objects rather than element coordinates, following function converts object element numbers to coordinates
    % Once the coordinates are determined, calculation of objects' dimensions and elimination of "surface objects" is much easier and faster
    [xx,yy,zz,brdr] = find_cond(cc.PixelIdxList{1,num},dims);

    % Dimensions (in voxels) of the considered object
    dim_X = max(xx)-min(xx)+1; 
    dim_Y = max(yy)-min(yy)+1; 
    dim_Z = max(zz)-min(zz)+1;

    % If the considered object touches one of the surfaces of the region of interest, or one of the object's dimension is greater or less than the maximum or minimum value, 
    % corresponding message is displayed and the program moves on to the next object. Otherwise, analysis is carried out and the object's geometric parameters are recorded.
    if brdr == 1
        fprintf(1,'%s%g%s\n','Pore # ',num,' is on the border ');   % Message: objects is on the surface of the ROI
        ii = ii + 1;                                                % Counter. Records how many "surface objects" have been detected
    elseif min([dim_X,dim_Y,dim_Z]) < pix_min
        fprintf(1,'%s%g%s%g%s\n','Pore # ',num,' has less than ', pix_min,...
            ' pixels across the smallest dimension');               % Message: object's smallest dimension is less than user specified value
        oo = oo + 1;                                                % Counter. Records how many "small objects" have been detected
    elseif max([dim_X,dim_Y,dim_Z]) > pix_max
        fprintf(1,'%s%g%s%g%s\n','Pore # ',num,' has more than ', pix_max,...
            ' pixels across the largest dimension');                % Message: object's largest dimension is greater than user specified value
        pp = pp + 1;                                                % Counter. Records how many "large objects" have been detected
    else

        % In this program 3D shapes are extracted using 'isosurface' function. To reduce processing time, the size of the input matrix for the isosurface function should be reduced.
        % This is done by creating a separate 3D matrix for the object. The dimensions of the matrix match the dimensions of the object.
        
        % Object is moved to the origin of coordinates so that it is contained within the new 3D matrix. The coordinates are shifted by '2' so that object starts at the second element
        % (there are no zeroth elements in Matlab. Starting with the first element will lead to object clipping by the 'isosurface' function)
        xx = xx - ones(size(xx))*min(xx) + 2;
        yy = yy - ones(size(yy))*min(yy) + 2;
        zz = zz - ones(size(zz))*min(zz) + 2;
        
        pore = zeros(dim_Y+2,dim_X+2,dim_Z+2);      % Memory allocation for the object's 3D matrix. The matrix's dimensions are equal to the object's dimension plus '2' (to prevent object clipping)

        % Object is copied to the new 3D matrix
        for i = 1:size(xx,1)
            pore(yy(i),xx(i),zz(i)) = 1;
        end

        poresurf_i = isosurface(pore,0.5);            % 3D shape is created based on the object's voxel coordinates
        
        % Number of faces in the extracted 3D shape
        facez = size(poresurf_i.faces,1);
        
        % The function 'isosurf_geom_params' computes the following geometric parameters:
        % - dimensions of the extracted 3D shape (dim_Xi, dim_Yi, dim_Zi)
        % - object's surface area (atotal), 
        % - volume (Vol), 
        % - center of mass (Xgtotal), 
        % - inertia tensor (J), 
        % - 3 principal moments of inertia (I11, I22, I33),
        % - semiaxes of the approximating ellipsoid (having the same principal moments of inertia and volume as the original object) - a,b,c
        % - axes of the principal moments of inertia (V)
        [dim_Xi,dim_Yi,dim_Zi,atotali,Voli,Xgtotali,Ji,I11i,I22i,I33i,ai,bi,ci,Vi] = isosurf_geom_params(poresurf_i);

        % The original object is: a) rotated so that principal axes are aligned with the coordinate system axes (I11 is aligned with x-axis, I22 - with y-axis, I33 - with z-axis)
        % b) resized so that object's largest dimension is equal to 2 (for FEA analysis)
        % c) moved so that object's midpoint lies at the origin of coordinates
        [poresurf_p,sf] = isosurf_transform(poresurf_i,Vi);
        
        % Geometric parameter computation for the transformed object
        [dim_Xp,dim_Yp,dim_Zp,atotalp,Volp,Xgtotalp,Jp,I11p,I22p,I33p,ap,bp,cp,Vp] = isosurf_geom_params(poresurf_p);

        fprintf(1,'%s%g%s%g%s%g%s%g%s%g%s%g%s%g\n','Pore # ',num,' is processed. Volume_i = ',Voli,'. I11i = ',I11i,'. I22i = ',I22i,'. I33i = ',I33i,'. Area_i = ',atotali,'. Faces: ',facez);

        % Select geometric parameters are stored in a summary table
        dat_table(m,1) = num; 
        dat_table(m,2:4) = [I11i, I22i, I33i];
        dat_table(m,5:7) = [ai, bi, ci];
        dat_table(m,8:10) = [Vi(1,1), Vi(2,1), Vi(3,1)];
        dat_table(m,11:13) = [Vi(1,2), Vi(2,2), Vi(3,2)];
        dat_table(m,14:16) = [Vi(1,3), Vi(2,3), Vi(3,3)];
        dat_table(m,17) = Voli;
        dat_table(m,18) = atotali;
        dat_table(m,19) = facez;
        dat_table(m,20) = sf;
        dat_table(m,21:23) = [dim_Xi, dim_Yi, dim_Zi];
        dat_table(m,24:26) = [dim_Xp, dim_Yp, dim_Zp];

        m = m+1;                                                        % Counter of sucessfully analyzed objects
        
        zeross = filename_app_zeros(5,num);                             % This function determines the number of zeros to be appended to the filename before the object number for consistent file numbering
        filename = strcat(init_filename, '_', zeross, num2str(num));    % Name of the object's files

        % The original and transformed 3D shapes are exported in STL format (triangular surface elements)
        stlwrite(strcat('stl\',filename,'_org.stl'),poresurf_i,'mode','binary')
        stlwrite(strcat('stl\',filename,'.stl'),poresurf_p,'mode','binary');

        % 4 views of the orginal 3D shape are plotted with global and principal axes and the image is then stored in PNG format
        set(0,'DefaultAxesPosition',[0.1 0.1 0.8 0.8],'DefaultFigureResize', 'on');
        figure
        subplot(2,2,1);
        p1 = patch(poresurf_i);
        set(p1,'facecolor',[0.85,1,0.85],'edgecolor','black');
        xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis'), lighting flat, material dull, camlight('headlight'), view([90 0]), axis image, grid on, box off
        axh = subplot(2,2,2);
        copyobj(p1,axh);
        xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis'), lighting flat, material dull, camlight('headlight'), view([0 0]), axis image, grid on, box off
        axh = subplot(2,2,3);
        copyobj(p1,axh);
        xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis'), lighting flat, material dull, camlight('headlight'), view([0 90]), axis image, grid on, box off
        axh = subplot(2,2,4);
        hold on, 
        copyobj(p1,axh);
        orgn = Xgtotali;
        V = Vi'; VV1 = V(1,1:3); VV2 = V(2,1:3); VV3 = V(3,1:3);
        % Global coordinate axes
        quiver3(orgn(1),orgn(2),orgn(3),1,0,0, dim_Xi*0.7, 'Color','r', 'LineWidth',1);
        quiver3(orgn(1),orgn(2),orgn(3),0,1,0, dim_Yi*0.7, 'Color','g', 'LineWidth',1);
        quiver3(orgn(1),orgn(2),orgn(3),0,0,1, dim_Zi*0.7, 'Color','b', 'LineWidth',1);
        % Principal axes
        quiver3(orgn(1),orgn(2),orgn(3),VV1(1),VV1(2),VV1(3), dim_Xp*sf*0.7, 'Color','r', 'LineWidth',1,'LineStyle','--');
        quiver3(orgn(1),orgn(2),orgn(3),VV2(1),VV2(2),VV2(3), dim_Yp*sf*0.7, 'Color','g', 'LineWidth',1,'LineStyle','--');
        quiver3(orgn(1),orgn(2),orgn(3),VV3(1),VV3(2),VV3(3), dim_Zp*sf*0.7, 'Color','b', 'LineWidth',1,'LineStyle','--');
        % Plot configuration
        xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis'), lighting flat, material dull, camlight('headlight'), view([60 20]), axis image, grid on, box off, hold off

        print('-dpng',strcat('png\',filename,'.png'))           % Export of the image to the PNG file

    end
    
    % The data table containing geometric parameters of the analyzed objects is stored in MAT format
    save(strcat(init_filename, '_processed', ext_mat),'dat_table')
end

% Summary of analysis
fprintf(1,'\n%s%g%s%g%s','- on the border: ',ii,' objects (',ii/ttl_num*100,'% of the total number) ');
fprintf(1,'\n%s%g%s%g%s%g%s','- have less than ', pix_min,' pixels across the smallest dimension: ',oo,' objects (',oo/ttl_num*100,'% of the total number) ');
fprintf(1,'\n%s%g%s%g%s%g%s','- have more than ', pix_max,' pixels across the largest dimension: ',pp,' objects (',pp/ttl_num*100,'% of the total number) ');
fprintf(1,'\n%s%g%s%g%s','Successfully written ',ttl_num-ii-oo-pp,' objects (',(ttl_num-ii-oo-pp)/ttl_num*100,'% of the total number). ');