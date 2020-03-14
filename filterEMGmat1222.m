function   EMG_filt_Mat = filterEMGmat1222(EMG_raw_Mat, hp1, lp1,lp2, timeBin)
%% fft of each filtering step
% hp1 = 35;
% lp1 = 50;
% lp2=6;
% timeBin = 10;

n = size(EMG_raw_Mat,2);    % number of columns
%% design filter
sf=1200;   %Hz - Sampling Freqeuncy
nf=sf/2;        %Hz - Nyquist Frequency - 1/2 Sampling Frequency
% design methods: passband>ellip
%HP1=designfilt('highpassiir','FilterOrder',5,'PassbandFrequency',hp1,...
%    'SampleRate',sf,'DesignMethod','butter');

% LP1=designfilt('lowpassiir','FilterOrder',5,'PassbandFrequency',lp1,...
%     'SampleRate',sf,'DesignMethod','butter');
% LP2=designfilt('lowpassiir','FilterOrder',4,'PassbandFrequency',lp2,...
%     'SampleRate',sf,'DesignMethod','butter');

HP1 = designfilt('highpassiir', 'FilterOrder', 5, 'PassbandFrequency', ...
    hp1, 'StopbandAttenuation', 60, 'PassbandRipple', 1, ...
    'SampleRate', 1200, 'DesignMethod', 'ellip');
if lp1~=0
    LP1 = designfilt('lowpassiir', 'FilterOrder', 5, 'PassbandFrequency', ...
        lp1, 'StopbandAttenuation', 60, 'PassbandRipple', 1, ...
        'SampleRate', 1200, 'DesignMethod', 'ellip');
end
LP2 = designfilt('lowpassiir', 'FilterOrder', 5, 'PassbandFrequency', ...
    lp2, 'StopbandAttenuation', 60, 'PassbandRipple', 1, ...
    'SampleRate', 1200, 'DesignMethod', 'ellip');


% fvtool(HP1)

%%
% EMG_raw_Mat = T1(3).EMG.Right(1).all;
n = size(EMG_raw_Mat,2);    % number of columns
% allMuscleNames = T1(3).EMG.Right(1).muscleName;

for i=1:n
    % 1)raw data
    temp1(:,i)= EMG_raw_Mat(:,i);
    % 2)High pass filter
    
    temp2(:,i)= filtfilt(HP1,temp1(:,i));
    if lp1~=0
        temp2(:,i)= filtfilt(LP1,temp2(:,i));
    end
    % 3)Detrend
    temp3(:,i)=detrend(temp2(:,i));
    
    % 4)Rectify
    temp4(:,i)=abs(temp3(:,i));
    
    % 5)Low pass filter
    
    temp5(:,i)= filtfilt(LP2,temp4(:,i));
    
    % 6)Devide by max value
    max_amp = max(temp5(:,i));
    %                max_amp=1;    % if division is not desired
    temp6(:,i) = temp5(:,i)/max_amp;
    % 7)time binning
    temp7(:,i) = (mean1(temp6(:,i)', timeBin,1))';
    
end
% Finally
EMG_filt_Mat = temp7;

% fft1(temp1, temp2, temp4, temp5)
%%
% for i = 1:n
%     handle = figure(i);
%     handle.Name = allMuscleNames{i};
%     fft1(temp1(:,i), temp2(:,i), temp4(:,i), temp6(:,i))
% %     suptitle(allMuscleNames{i})
%     
% end

% for i = 1:n
%     figure(i)
%     suptitle(allMuscleNames{i})
%     
% end

%% fft function for plot
% function fft1(raw, hp, rec, lp)
% allTitles = {'Raw' 'HP+LP' 'Rectified' 'LP'};
% % figure
% Fs = 1200;            % Sampling frequency
% TT = 1/Fs;             % Sampling period
% 
% k=1;
% for S = {raw, hp, rec, lp}
%     
%     S = cell2mat(S);
%     L = length(S);             % Length of signal
%     t = (0:L-1)*TT;        % Time vector
%     Y = fft(S);
%     P2 = abs(Y/L);
%     P1 = P2(1:L/2+1);
%     P1(2:end-1) = 2*P1(2:end-1);
%     
%     f = Fs*(0:(L/2))/L;
%     
%     subplot(4,2,2*k-1); hold on;
%     plot(f,P1)
%     title(allTitles{k})
%     
%     subplot(4,2,2*k);   hold on;
%     plot(S/max(S))
%     
%     k=k+1;
% end
% 
% 
% 
% 
% % suptitle('Single-Sided Amplitude Spectrum of S(t)')
% 
% 
% 
% 
% end
% 


