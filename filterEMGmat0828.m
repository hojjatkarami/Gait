function   EMG_filt_Mat = filterEMGmat0828(EMG_raw_Mat, highPassFreq, lowPassFreq, timeBin)

% each column contains a signal to be filtered
n = size(EMG_raw_Mat,2);    % number of columns
%% design filter
sf=1200;   %Hz - Sampling Freqeuncy
fn=sf/2;        %Hz - Nyquist Frequency - 1/2 Sampling Frequency
bin = timeBin;

HP30=designfilt('highpassiir','FilterOrder',5,'HalfPowerFrequency',highPassFreq,'SampleRate',1200,'DesignMethod','butter');
LP10=designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',lowPassFreq,'SampleRate',1200,'DesignMethod','butter');
%%

for i=1:n
    % 1)raw data
    temp1(:,i)= EMG_raw_Mat(:,i);
    % 2)High pass filter
    temp2(:,i)= filtfilt(HP30,temp1(:,i));
% temp2(:,i)=temp1(:,i);
    % 3)Detrend
    temp3(:,i)=detrend(temp2(:,i));
    % 4)Rectify
    temp4(:,i)=abs(temp3(:,i));
    % 5)Low pass filter
    temp5(:,i)= filtfilt(LP10,temp4(:,i));
    % 6)Devide by max value
    max_amp = max(temp5(:,i));
%                max_amp=1;    % if division is not desired
    temp6(:,i) = temp5(:,i)/max_amp;
    % 7)time binning
    temp7(:,i) = (mean1(temp6(:,i)', bin,1))';
end
% Finally
EMG_filt_Mat = temp7;  
    

