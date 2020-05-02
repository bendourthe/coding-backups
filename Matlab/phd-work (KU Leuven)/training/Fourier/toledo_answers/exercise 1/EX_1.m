%% 2. EMG signaal vermoeidheid verwerken
%---------------------------------------
clear all; close all; clc;
load('Exercise1.mat');

%% 1. plot the emg signal
figure()
plot(Delt_l,'b')
title('activity Deltoideus')
xlabel('frames')
ylabel('mV')

%% 2. offset correction
Delt_l=Delt_l-mean(Delt_l);

%% 3. plot the frequency spectrum of the left Delt during the whole movement
Delt_l_ft=fft(Delt_l);
sze = length(Delt_l);
ff= fix(sze/2) + 1; 
f = [0:ff-1]*fs/sze;
figure
plot(f(1:ff), abs(Delt_l_ft(1:ff)/sze*2));
xlabel('Frequency in Hz');
ylabel('Magnitude');
title('Frequency spectrum');
axis tight;

% je kan dit ook doen met de functie die gegeven is (make_freq_spect)
make_freq_spect(Delt_l,fs,1);


%% 4. Bandstop filter

cutoff=[49 51];
order=4;
[b, a] = butter(order, cutoff/(0.5*fs),'stop');
Delt_l_filt = filtfilt(b, a, Delt_l); 

% plot the frequency spectrum
[~, bandx_data, bandy_data]=make_freq_spect(Delt_l_filt,fs,1);
title(['Frequency spectrum after bandpass filtering ' num2str(cutoff)])

%% 5. Calculate the average frequency and activity for each second

Gemfreq=zeros(fix(length(Delt_l_filt)./fs),1);% pre allocatie van de variabelen
Mean_amp=zeros(fix(length(Delt_l_filt)./fs),1);% pre allocatie van de variabelen

% opmerkinge: het pre alloceren van de variablen is niet essentieel, maar
% dit versnelt de lus wel. Indien je de variabelen nog niet aanmaakt,
% veranderd de lengte van de variabele in elke lus. 

        %% 5.1: Eerste manier om lus te vormen

% Calculate the mean activity and mean frequency for each second
% We maken dus een lus die in stappen van 1 gaan en start bij 1. De lus
% moet lopen tot de laatste volledige seconde. Het aantal seconden van het
% signaal is gelijk aan het aantal frames (length(Delt_l_filt)) gedeeld
% door de sampling frequency (fs). We ronden dit af naar beneden om
% uitsluitend volledige seconden te hebben. (fix gebruiken i.p.v. round)

for i=1:fix(length(Delt_l_filt)./fs)
    % selecteer de frames (vb 1 tot 2000, 2001 tot 4000, enz...)
    time_sel=((i-1)*fs)+1:(i)*fs;    
    
    % selecteer de data voor de geselecteerde frames
    data=Delt_l_filt(time_sel);
    
    % we willen geen grafieken plotten omdat dit teveel zou vergen van het
    % geheugen. Daarom maken we de bool_plot=0.
    bool_plot=0;
   
    % get the frequency spectrum
    [data_ft,x_data,y_data] = make_freq_spect( data,fs,bool_plot);
    
    % bereken de gemiddelde frequentie. Dit is een gewogen gemiddelde. De
    % x_data zijn alle frequenties (frequentie van de harmischen). De y-data
    % is de magnitude(power: amplitude harmonischen). Met andere woorden de
    % y-data bepaalt hoe sterk die frequentie doorweegt in het signaal.
    % Bijgevolg wordt de som genomen van de y_data*x_data. Deze som wordt
    % dan gedeeld door de sum van de y_data. 
    
    Gemfreq(i) = sum(y_data.*x_data)./sum(y_data);
    % get the mean amplitude
    
    % dit berekenen we door het gemiddelde te nemen van de absolute waarde.
    % We nemen de absolute waarde omdat we het signaal willen
    % gelijkrichten. 
    % Het nemen van de absolute waarde van een emg signaal word 
    % ook wel rectificatie genoemd.
    
    Mean_amp(i)=mean(abs(data));    
end

figure()
subplot(2,1,1)
plot(Gemfreq,'r*');
xlabel('time(s)');
ylabel('Mean Frequency');
title('The influence of fatigue on mean frequency')

subplot(2,1,2)
title('The influence of fatigue on  mean amplitude')
plot(Mean_amp,'r*');
xlabel('time (s)');
ylabel('Mean amplitude');

        %% 5.2 tweede manier om lus te vormen
       
        
% op deze manier maken we een lus die in stappen van de de sampling frequency gaat
% in dit geval in stappen van 2000. Dus lus moet lopen tot het laatste
% veelvoud van 2000 dat kleiner of gelijk is aan lengte van het signaal. Bijgevolg
% nemen we de lengte van het signaal, en berekenden we de rest van de
% deling van de lengte door de sampling frequency. Dit geeft ons het
% laatste mogelijke veelvoud van 2000.

count=1;

for i=1:fs:(length(Delt_l_filt)-mod(length(Delt_l_filt),fs))
    time_sel= i:(i+2000);    
    data=Delt_l_filt(time_sel);
    bool_plot=0;
    % get the frequency spectrum
    [data_ft,x_data,y_data] = make_freq_spect( data,fs,bool_plot);
    % get the mean frequency
    Gemfreq(count) = sum(y_data.*x_data)./sum(y_data);
    % get the mean amplitude
    Mean_amp(count)=mean(abs(data));  
    count=count+1;
end

% hierbij maken we dus gebruik van een teller (count).
