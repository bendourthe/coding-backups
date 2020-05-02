clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DynamicCT_Prox: use the ProxFunction on multiple stl files to calculate:
%   - minimal joint space for each frame
%   - proximity areas for each joint
%   - displacement of the centroid of each joint proximity area
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input: 
%           - stl files of each bone
%           - number of frames
%           - proximity parameters
%           - azimuth and elevation for 3D visualization
% 
% Output:
%           - saves each figure in the current directory
%           - saves the minimal joint space and proximity areas in text
%           files:
%               'TMC Joint': minimal joint space (mm), proximity
%               area MC1 (mm2), proximity area trapezium (mm2)
%               'STT Joint': minimal joint space (mm), proximity
%               area scaphoid (mm2), proximity area trapezium (mm2)
%
% Dependencies:        
%               STL_ReadFile.m 
%                   TRI_RemoveInvalidTriangles.m
%                   TRI_RemoveBadlyConnectedTriangles.m
%               GUI_PlotShells
%               ProxFunction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Selection of the current directory (where the STL files are)
cd 'H:\Data Faes\Patient_07_Post\STLs\4D_AutoSeg_Output_Flex\stls\'
dir = 'H:\Data Faes\Patient_07_Post\STLs\4D_AutoSeg_Output_Flex\stls\';

% Name of the text files where the numerical data will be saved
filenameTMC = 'TMC_joint.txt';
    fileTMC = [dir filenameTMC];
filenameSTT = 'STT_joint.txt';
    fileSTT = [dir filenameSTT];
    
% Number of frames
num_frames = 25;

% Maximal distance (defines the maximal range under which the calculation
% of distance is done)
d_max = 3;       % in mm
% Proximity distance (defines the minimal joint space under which we can
% expect deformation of the cartilage)
d_prox = 1.5;    % in mm

% Orientation of the figures (define azimith and elevation for each plot)
    % MC1
az1 = 230;
el1 = -30;
    % Trap (TMC)
az2 = -3200;
el2 = 60;
    % Scap
az3 = -280;
el3 = 60;
    % Trap (STT)
az4 = 230;
el4 = -30;

% FILES
    MC1 = num2str([1:num_frames].','Patient_07_4D_Flex_Frame25_MC1_SE000000_Frame%04d.stl');
    Trap = num2str([1:num_frames].','Patient_07_4D_Flex_Frame25_Trap_SE000000_Frame%04d.stl');
    Scap = num2str([1:num_frames].','Patient_07_4D_Flex_Frame25_Scap_SE000000_Frame%04d.stl');

% TMC Joint:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:num_frames;
    [TMC_dm(i),TMC_PA_mc1(i),TMC_PA_trap(i),TMC_ctr1(i),TMC_ctr2(i),TMC_ctr3(i),area1,area2] = ...
        ProxFunction(dir,d_max,d_prox,MC1(i,:),Trap(i,:),az1,el1,az2,el2);
    set(area1,'numbertitle','off','name',sprintf('TMC joint: Proximity pattern MC1 #%d',i));
    % Save proximity map of the MC1
    savefig(area1,sprintf('TMC_Prox_MC1_#%d.fig',i));
    set(area2,'numbertitle','off','name',sprintf('TMC joint: Proximity pattern Trapezium #%d',i));
    % Save proximity map of the Trapezium
    savefig(area2,sprintf('TMC_Prox_Trap_#%d.fig',i));
end

% STT Joint:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:num_frames;
    [STT_dm(i),STT_PA_scap(i),STT_PA_trap(i),STT_ctr1(i),STT_ctr2(i),STT_ctr3(i),area1,area2] = ...
        ProxFunction(dir,d_max,d_prox,Scap(i,:),Trap(i,:),az3,el3,az4,el4);
    set(area1,'numbertitle','off','name',sprintf('STT joint: Proximity pattern Scaphoid #%d',i));
    % Save proximity map of the Scaphoid
    savefig(area1,sprintf('STT_Prox_Scap_#%d.fig',i));
    set(area2,'numbertitle','off','name',sprintf('STT joint: Proximity pattern Trapezium #%d',i));
    % Save proximity map of the Trapezium
    savefig(area2,sprintf('STT_Prox_Trap_#%d.fig',i));
end

% Save minimal distances and proximity areas (mm2) in text files
    % TMC Joint    
    fid = fopen(fileTMC,'w+');
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\r\n',[TMC_dm;TMC_PA_mc1;TMC_PA_trap]);
    fclose(fid);
    % STT Joint
    fid = fopen(fileSTT,'w+');
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\r\n',[STT_dm;STT_PA_scap;STT_PA_trap]);
    fclose(fid);
    
% Plots the Trapezium (frame 1) with the displacement of the centroid of
% the proximity area

% Read STL files
    [F_Trap, V_Trap] =  STL_ReadFile([dir Trap(1,:)],true);

% TMC Joint    
    f1 = figure;
    [obj, li, ax] = GUI_PlotShells(f1, {F_Trap},{V_Trap},...
        {ones(size(V_Trap,1),1)});
    box off
    view(az2,el2);
    hold on
    for i=1:num_frames
        plot3(TMC_ctr1(i),TMC_ctr2(i),TMC_ctr3(i),'Marker','o', ...
            'MarkerSize',10,'Color','k','MarkerFaceColor','k');
        text(TMC_ctr1(i),TMC_ctr2(i),TMC_ctr3(i),sprintf('#%d',i),'FontSize',15,'FontWeight','bold');
    end
    set(f1,'numbertitle','off','name','TMC Joint: Proximity centroid displacement');
    % Save figure
    savefig(f1,'TMC_PA_ctr_mot');
    
% STT Joint
    f2 = figure;
    [obj, li, ax] = GUI_PlotShells(f2, {F_Trap},{V_Trap},...
        {ones(size(V_Trap,1),1)});
    box off
    view(az4,el4);
    hold on
    for i=1:num_frames
        plot3(STT_ctr1(i),STT_ctr2(i),STT_ctr3(i),'Marker','o', ...
            'MarkerSize',10,'Color','k','MarkerFaceColor','k');
        text(STT_ctr1(i),STT_ctr2(i),STT_ctr3(i),sprintf('#%d',i),'FontSize',15,'FontWeight','bold');
    end
    set(f2,'numbertitle','off','name','STT Joint: Proximity centroid displacement');
    % Save figure
    savefig(f2,'STT_PA_ctr_mot');