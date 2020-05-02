clear all; close all;
dir = 'C:\Users\u0096636\Documents\Professional\PhD KU Leuven\Data\Forward Model\Results\New Implementation\';
filenameRot = '0N_rot_acc.fig';
    fileRot = [dir filenameRot];
open([dir filenameRot]);
D=get(gca,'Children'); %get the handle of the line object
XData=get(D,'XData'); %get the x data
YData=get(D,'YData'); %get the y data

t = cell2mat(XData(1,1));
Rot_x = cell2mat(YData(4,1));
Rot_y = cell2mat(YData(3,1));
Rot_z = cell2mat(YData(2,1));

for i=1:length(Rot_x)
    Rot_acc_x(i) = sum(Rot_x(:, 1:i));
    Rot_acc_y(i) = sum(Rot_y(:, 1:i));
    Rot_acc_z(i) = sum(Rot_z(:, 1:i));
end

f1=figure;
plot(t,Rot_x,'r',t,Rot_y,'g',t,Rot_z,'b',t,t*0,'k');
xlabel('Time (s)');
ylabel('Rotation angle (degrees)')
set(gcf,'numbertitle','off','name','Evolution of the rotation angles');

f2=figure;
plot(t,Rot_acc_x,'r',t,Rot_acc_y,'g',t,Rot_acc_z,'b',t,t*0,'k');
xlabel('Time (s)');
ylabel('Rotation angle accummulated (degrees)')
legend('Rot_x','Rot_y','Rot_z','Location','NorthEastOutside')
set(gcf,'numbertitle','off','name','Evolution of the rotation angles (accumulated)');