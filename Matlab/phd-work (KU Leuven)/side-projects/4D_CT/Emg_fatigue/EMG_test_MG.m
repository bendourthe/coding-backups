clear all; close all;clc;
%%
DATA = loadmat('Fatigue2682013.sig');
Delt_l = DATA(:,1);
Obliq_r = DATA(:,2);
Delt_r = DATA(:,3);


for j = 1:3;
    % Offset correction
    data(:,j) = DATA(:,j) - mean(DATA(:,j));
    
    for i = 1:4;
        fs = 2000; % sample ferquency
    
        % Divide in 4 windows (5000Hz).
        data_div(:,i) = data(((i*5000)-5000)+1:(i*5000));

        % Determine frequency distribution.
        [Pxx(:,i),Fx(:,i)] = pwelch(data_div(:,i),1024,512,500,fs);

        % Determine medianfrequency
        Medfreq(:,i) = Fx(sum(cumsum(Pxx(:,i))<= sum(Pxx(:,i))./2));

        % Determine meanfrequency.
        Gemfreq(:,i) = sum(Pxx(:,i).*Fx(:,i))./sum(Pxx(:,i));

        % Mean amplitude after rectifying EMG.
        data_abs(:,i) = abs(data_div(:,i));
        data_mean(:,i) = mean(data_abs(:,i));

        figure(j)
        plot(Fx(:,j), Pxx(:,j));
        title 'Powerspectrum '
        xlabel 'Frequency [Hz]'
        ylabel 'Power'
    end

        figure(j+3);
        subplot(3,1,1)
        hold on
        bar(Medfreq);
        title 'Median frequency for the 4 windows'
        xlabel 'Windows'
        ylabel 'Median frequency [Hz]'

        subplot(3,1,2);
        bar(Gemfreq, 'r');
        title 'Mean frequency for the 4 windows'
        xlabel 'Windows'
        ylabel 'Mean frequency [Hz]'

        subplot(3,1,3);
        bar(data_mean, 'g');
        title 'Mean amplitude'
        xlabel 'Windows'
        ylabel 'Mean amplitude [Hz]'

       

    %     wm=400/2000/2;
    %     [a,n]=butter(2,wm,'low');
    %     filtdata2(:,i)=filtfilt(a,n,filtdata(:,i));

    %     figure()
    %     fs = 2000;
    %     t = 0:1/fs:2;
    %     x = filtdata;
    %     spectrogram(x,kaiser(256,5),220,512,fs,'yaxis')
    %     ylim([0 250]);

%         figure(2)
%         test2=pwelch(filtdata(1:20000,i));
%         subplot(3,1,i), plot(test2)
%         xlim([0 750]);
%         hold on
%         second2=pwelch(filtdata(end-20000:end,i));
%         plot(second,'r')

end