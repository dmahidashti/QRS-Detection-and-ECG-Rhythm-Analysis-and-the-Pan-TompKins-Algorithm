%ECG Beat Detection and the Pan-TompKins Algorithm

%% Load Data

ECG3 = load('/Volumes/GoogleDrive/My Drive/University/BME772 Signal Analysis/Lab 3/Data/ECG3.txt');
ECG4 = load('/Volumes/GoogleDrive/My Drive/University/BME772 Signal Analysis/Lab 3/Data/ECG4.txt');
ECG5 = load('/Volumes/GoogleDrive/My Drive/University/BME772 Signal Analysis/Lab 3/Data/ECG5.txt');
ECG6 = load('/Volumes/GoogleDrive/My Drive/University/BME772 Signal Analysis/Lab 3/Data/ECG6.txt');

%% Graphing Original ECG Signals
fs = 200;
%Time Vectors
t = (0:length(ECG3)-1)/fs;

%Plots
figure;
subplot(4,1,1);
plot(t, ECG3);
title('Original ECG3');
xlabel('Time (sec)');
ylabel('Amplitude (mV)');
legend('ECG Signal');
subplot(4,1,2);
plot(t, ECG4);
title('Original ECG4');
xlabel('Time (sec)');
ylabel('Amplitude (mV)');
legend('ECG Signal');
subplot(4,1,3);
plot(t, ECG5);
title('Original ECG5');
xlabel('Time (sec)');
ylabel('Amplitude (mV)');
legend('ECG Signal');
subplot(4,1,4);
plot(t, ECG6);
title('Original ECG6');
xlabel('Time (sec)');
ylabel('Amplitude (mV)');
legend('ECG Signal');

%% Defining Filter Coefficients, Fequency Response and ZPlanes
M=128;
%Low Pass Filter
lp_b = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
lp_a = [1 -2 1 0 0 0 0 0 0 0 0 0 0];
figure;
freqz(lp_b,lp_a,M,fs);
title('Frequency Response of Low Pass Filter');
figure;
zplane(lp_b,lp_a);
title('Location of Poles and Zeros of Low Pass Filter');

%High Pass Filter
hp_b = [-1/32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1/32];
hp_a = [1 -1];
figure;
freqz(hp_b,hp_a,M,fs);
title('Frequency Response of High Pass Filter');
figure;
zplane(hp_b,hp_a);
title('Location of Poles and Zeros of High Pass Filter');

% Bandpass Filter
bp_b=conv(lp_b,hp_b);
bp_a=conv(lp_a,hp_a);
figure;
freqz(bp_b,bp_a,M,fs);
title('Frequency Response of Bandpass Filter');
figure;
zplane(bp_b,bp_a);
title('Location of Poles and Zeros of Bandpass Filter');

%Derivative Based Filter
db_b = [0.25 0.125 0 -0.125 -0.25];
db_a = [1 0 0 0 0];
figure;
freqz(db_b,db_a,M,fs);
title('Frequency Response of Derivative Based Filter');
figure;
zplane(db_b,db_a);
title('Location of Poles and Zeros of Derivative Based Filter');

%Moving Window Filter
mw_b = [1:30];
mw_b(1:30)=1;
mw_b = (1/30)*mw_b;
mw_a = 1;
figure;
freqz(mw_b,mw_a,M,fs);
title('Frequency Response of Moving Window Filter');
figure;
zplane(mw_b,mw_a);
title('Location of Poles and Zeros of Moving Window Filter');

%% Specify the Input ECG -> ECG3, ECG4, ECG5, ECG6

ECG = ECG3; 

M1= ECG;
E1=['Error'];
if M1==ECG3
    E1=['ECG3'];
else if M1==ECG4
        E1=['ECG4'];
    else if M1==ECG5
            E1=['ECG5'];
        else if M1==ECG6
                E1=['ECG6'];
            end
        end 
    end
end

t = (0:length(ECG)-1)/fs;
%% Building Filters and filtering the signal + squaring operation

% Bandpass Filter
bp_ecg = filter(bp_b,bp_a,ECG);
% Derivative Based Filter
db_ecg = filter(db_b,db_a,bp_ecg);
% Squaring Operation
sq_ecg = db_ecg.^2;
% Moving Window Filter
mw_ecg = filter(mw_b,mw_a,sq_ecg);

%Plotting
figure;
subplot(4,1,1);
plot(t,bp_ecg);
title(E1,'Signal After Bandpass Filter');
xlabel('Time (sec)');
ylabel('Amplitude (mV)');
subplot(4,1,2);
plot(t,db_ecg);
title('Signal After Derivative based Filter');
xlabel('Time (sec)');
ylabel('Amplitude (mV)');
subplot(4,1,3);
plot(t,sq_ecg);
title('Signal After Squaring Operation');
xlabel('Time (sec)');
ylabel('Amplitude (mV)');
subplot(4,1,4);
plot(t,mw_ecg);
title('Signal After Moving Window Filter');
xlabel('Time (sec)');
ylabel('Amplitude (mV)');

%% Peak Detection
ECG_filtered = mw_ecg;

% Starting Point is at time=0.37 which is at index=74, re-align to start at index 74
ECG_aligned = mw_ecg(74:length(mw_ecg));
time_aligned = t(74:length(t));

% Use Local Maxima Finder to locate the peaks
[pks,locs] = findpeaks(ECG_aligned,time_aligned,'MinPeakDistance',0.4,'MinPeakHeight',1000000);
 
% Plotting QRS Complex Peaks
figure;
hold on
findpeaks(ECG_aligned,time_aligned,'MinPeakDistance',0.4,'MinPeakHeight',1000000);
hold off
title(E1,'Peaks of the QRS Complex for the ECG Signal');
xlabel('Time (sec)');
ylabel('Amplitude (mV)');
legend('ECG Signal','QRS Complex Peak');

%% RR Interval Calculation

%Finding the difference between the peak locations
rr_interval=diff(locs);
mean_rr_interval=mean(rr_interval)*1000;
std_rr_interval=std(rr_interval)*1000;
rr_disp =['In ',E1,' the mean of RR intervals is ',num2str(mean_rr_interval),' (ms) and the standard deviation is ',num2str(std_rr_interval),' (ms)'];
disp(rr_disp);

%% Heart Rate Calculation

peak_count=length(pks);
%peak_time = length(time_aligned)/fs;
%bpm =(peak_count/peak_time)*60;
bpm=60./(mean(rr_interval));
bpm_disp =['In ',E1,' there are ',num2str(peak_count),' beats in total, the beats per minute is ',num2str(bpm)];
disp(bpm_disp);

%% QRS Complex Width

% Setup to cater to the specific signal
if M1==ECG3
   ECG_filtered=ECG_filtered(65:length(ECG_filtered)-1);
   distance=50;
else if M1==ECG4
       ECG_filtered=ECG_filtered(65:length(ECG_filtered)-1);
       distance=40;
    else if M1==ECG5
           ECG_filtered=ECG_filtered(200:length(ECG_filtered)-1);
           distance=40;
        else if M1==ECG6
                ECG_filtered=ECG_filtered(200:length(ECG_filtered)-1);
                distance=90;
            end
        end 
    end
end 

% Redefining time
t_new = (0:length(ECG_filtered)-1)/fs;
% Finding the local minima => Q and S points
[TF,P] = islocalmin(ECG_filtered,'MinSeparation',t_new(distance),'SamplePoints',t_new);
% Plotting positions of the Q and S point
figure;
%subplot(2,1,1);
hold on
findpeaks(ECG_filtered,t_new,'MinPeakDistance',0.4,'MinPeakHeight',1000000);
plot(t_new(TF), ECG_filtered(TF),'r*')
hold off
title(E1,'Filtered ECG with Q and S Points');
xlabel('Time (sec)');
ylabel('Amplitude (mV)');
legend('ECG Signal','QRS Complex Peak','Q and S Points of the Complex');

%Finding the QRS Width
Q=zeros(1,60);
i=1;
while i<60
    Q(i)= (TF(i+2)-TF(i))/fs;
    i = i+2;
end
% Averaging the Widths
Avg_QS=mean(Q);
% Displaying the average QRS width
QS_disp =['In ',E1,' the mean QRS Complex Width is ',num2str(Avg_QS)];
disp(QS_disp);



