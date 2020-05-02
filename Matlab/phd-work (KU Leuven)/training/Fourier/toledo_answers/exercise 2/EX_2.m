%% 1. load the data
load('Exercise2.mat');
fs=1000;
figure()
% plot het signaal
plot(Gastrocn);
title('activity of the Gastrocn during walking')

% plot het frequentie spectrum
[Gastrocn_ft,rawx_data,rawy_data] = make_freq_spect(Gastrocn,fs,1);

%% 2. implement a bandpass filter

cutoff=[50 400];
order=4;
fs=1000;
[b, a] = butter(order, cutoff/(0.5*fs));
Gastrocn_band = filtfilt(b, a, Gastrocn); 
[~, bandx_data, bandy_data]=make_freq_spect(Gastrocn_band,fs,0);
title(['Frequency spectrum after bandpass filtering ' num2str(cutoff)])

% compare both frequency spectra
figure()
subplot(1,2,1)
title('frequency spectrum raw signal')
plot(rawx_data,rawy_data);
subplot(1,2,2)
plot(bandx_data,bandy_data);
title('frequency spectrum band pass filtered signal')

%% 3. implement a lowpass filter

cutoff=200;
order=4;
fs=1000;
[b, a] = butter(order, cutoff/(0.5*fs));
Gastrocn_low = filtfilt(b, a, Gastrocn); 
[~, lowx_data, lowy_data]=make_freq_spect(Gastrocn_low,fs,0);
title(['Frequency spectrum after bandpass filtering ' num2str(cutoff)])

% compare both frequency spectra
figure()
subplot(1,2,1)
title('frequency spectrum raw signal')
plot(rawx_data,rawy_data);
subplot(1,2,2)
plot(lowx_data,lowy_data);
title('frequency spectrum low pass filtered signal')

%% 4. implement a highpass filter

cutoff=100;
order=4;
fs=1000;
[b, a] = butter(order, cutoff/(0.5*fs),'high');
Gastrocn_high = filtfilt(b, a, Gastrocn); 
[temp, highx_data, highy_data]=make_freq_spect(Gastrocn_high,fs,0);
title(['Frequency spectrum after bandpass filtering ' num2str(cutoff)])

% compare both frequency spectra
figure()
subplot(1,2,1)
title('frequency spectrum raw signal')
plot(rawx_data,rawy_data);
subplot(1,2,2)
plot(highx_data,highy_data);
title('frequency spectrum high pass filtered signal')

%% 5. Compare 