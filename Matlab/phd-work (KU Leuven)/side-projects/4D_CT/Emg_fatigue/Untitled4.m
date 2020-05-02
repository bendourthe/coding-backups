%
data=loadmat('Fatigue2682013.sig');

gem_freq1=average_freq_emg(data(:,1),121);
gem_freq2=average_freq_emg(data(:,2),121);
gem_freq3=average_freq_emg(data(:,3),121);

figure
plot(gem_freq1','*r')
hold on
plot(gem_freq2','*b')
plot(gem_freq3','*g')
xlabel('time (s)')
ylabel('mean frequency')
legend('R m. Deltoideus','R m.Obliquus','L m. Deltoideus')
