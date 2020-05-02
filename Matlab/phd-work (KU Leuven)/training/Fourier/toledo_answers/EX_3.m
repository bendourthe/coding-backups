%% Filter the gong
%-----------------
clear all; close all; clc
load gong.mat;
gong=y;
disp('GONG!')
[ft_out,x,y]=make_freq_spect(gong,Fs,1);
sound(gong,Fs)

% apply high pass filter on the signal
order=4;
cutoff=1000;
[b, a] = butter(order, cutoff/(0.5*Fs),'high');
gong_high_pass = filtfilt(b, a,gong);  
[ft_out,x,y]=make_freq_spect(gong_high_pass,Fs,1);
disp('GONG high pass filtered!')
sound(gong_high_pass,Fs);

% apply low pass filter on the signal
order=4;
cutoff=1000;
[b, a] = butter(order, cutoff/(0.5*Fs),'low');
gong_low_pass = filtfilt(b,a,gong);  
[ft_out,x,y]=make_freq_spect(gong_low_pass,Fs,1);
disp('GONG low pass filtered!')
sound(gong_low_pass,Fs);