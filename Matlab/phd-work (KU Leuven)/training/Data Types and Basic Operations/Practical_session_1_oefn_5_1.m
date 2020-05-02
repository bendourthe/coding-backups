% EMG processing

%1. select sit to stand and stand to sit
Sit_to_stand_time=Vastus_lateralis_raw(1500:3000,1);
Sit_to_stand_raw=Vastus_lateralis_raw(1500:3000,2);

Stand_to_sit_time=Vastus_lateralis_raw(4700:6100,1);
Stand_to_sit_raw=Vastus_lateralis_raw(4700:6100,2);

% opmerking het is ook mogelijk om zowel de tijd als de waarden in één
% matrix te laten staan. Maarde tijd en waarde splitsen in 2 vectoren is
% misschien duidelijker, maar wel meer typ werk. Beslis dus zelf hoe je dit
% probleem wilt aanpakken. 

% m.b.t. het selecteren van de datapunten om de beweging te starten en te
% stoppen: Nu is dit een zoektocht naar de juist indices, later zullen we
% zien hoe we dit met logisch indexeren sneller en correcter kunnen doen
% (de functie find of > <).

%2. rectification

Sit_to_stand_rect=abs(Sit_to_stand_raw);
Stand_to_sit_rect=abs(Stand_to_sit_raw);

%3. calculate mean

mean_activity_to_stand=mean(Sit_to_stand_rect);
mean_activity_to_sit=mean(Stand_to_sit_rect);

% 4. normalize to max value

max_to_stand=max(Sit_to_stand_rect);
max_to_sit=max(Stand_to_sit_rect);
max_both=max([max_to_stand max_to_sit]); % get the max value of both movement

Sit_to_stand_rect_norm=Sit_to_stand_rect./max_both;
Stand_to_sit_rect_norm=Stand_to_sit_rect./max_both;

%5. time maximally activated

[max_to_stand location_max_stand]=max(Sit_to_stand_rect); % search for the index (location_max_stand)
time_max_stand=Sit_to_stand_time(location_max_stand);     % search for the time frame based on the calculated index.

[max_to_sit location_max_sit]=max(Stand_to_sit_rect);
time_max_sit=Stand_to_sit_time(location_max_sit);

% 6. peak values + percentage movement

% max values already calculated ( max_to_stand max_to_sit)
% time values already calculated (time_max_stand time_max_sit)

% percentage= time from start to peak value/ time movement
time_start_stand_to_peak=time_max_stand-time_sit_to_stand(1);
perc_peak_stand=time_start_stand_to_peak./(time_sit_to_stand(2)-time_sit_to_stand(1));


time_start_sit_to_peak=time_max_sit-time_stand_to_sit(1);
perc_peak_sit=time_start_sit_to_peak./(time_stand_to_sit(2)-time_stand_to_sit(1));







