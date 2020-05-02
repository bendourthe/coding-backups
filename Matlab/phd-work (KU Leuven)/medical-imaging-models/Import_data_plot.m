clear all; close all;
dir = 'C:\Users\u0096636\Documents\Professional\PhD KU Leuven\Data\Forward Model\Results\New Implementation\';
filenameRot = '0N_Tz.fig';
    fileRot = [dir filenameRot];
open([dir filenameRot]);
D=get(gca,'Children'); %get the handle of the line object
XData=get(D,'XData'); %get the x data
YData=get(D,'YData'); %get the y data

t = cell2mat(XData(1,1));
Tin_x = cell2mat(YData(4,1));
Tlig_x = cell2mat(YData(3,1));
Tcart_x = cell2mat(YData(2,1));

Tcart_x = Tcart_x + Tin_x;
Tin_x = Tin_x - Tin_x;

f1=figure;
plot(t,Tin_x,'r',t,Tlig_x,'g',t,Tcart_x,'b');
xlabel('Time (s)');
ylabel('Tz (N.mm)');
legend('Tin-z','Tlig-z','Tcart-z','Location','NorthEastOutside');
set(gcf,'numbertitle','off','name','Evolution of the torque on Z');