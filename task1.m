
%% creating signals a)
fs=44100;
l=2;

%sample number
n=1:l*fs;     %unit sample 

%time vector
t_time=1/fs*(n-1);    %unit s



%% creating signals b)

y_unit=[ones(1,1),zeros(1,l*fs-1)];
stem(t_time,y_unit);
title('Impulse signal');
xlabel('time/s');
ylabel('amp');


%% creating signals c)
f1=1000;
f2=200;
a1=0.5;
a2=0.5;
t_s=1/(2*fs)*n(1:2*fs);
t_s=0:1/fs:2;
x1=a1*cos(2*pi*f1*t_s);
x2=a2*cos(2*pi*f2*t_s);
xsum=x1+x2;
figure;
plot(t_s,xsum);
title('The sum of two cos signals');
xlabel('sample num');
ylabel('amp');
sound(xsum,fs);    %can hear two impulse sound

%% creating signals d)
[audio,Fs]=audioread("audio\audio\drumloop.wav");
fprintf('the default sampeling rate is:%d\n',Fs);

%=======hear the audio========%
sound(audio,Fs);
figure;
plot(t_time,audio(n));
title('The amplitude of drum signal');
xlabel('time/s');
ylabel('amp');

%% creating signals e)
b = [1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8];
a=1;

y_filter = filter(b,a,audio);

figure;
plot(t_time,audio(n));
hold on;
plot(t_time,filter(b,a,audio(n)));
legend('Input Data','Filtered Data');
title('Audio signal after Using FIR filter comparision');
xlabel('time/s');
ylabel('amp');
%write the filter audiio file into wav file
audiowrite('task1_audio.wav',y_filter,Fs);
sound(y_filter,Fs);
%This is a low pass filter
%% creating signals f)
figure;
y_response = filter(b,a,y_unit);
plot(t_time,y_response);
title('Impulse response signal');
xlabel('time/s');
ylabel('amp');

%% Plotting a)
figure;
stem(t_time,y_unit);
hold on;
stem(t_time,y_response);
legend('Impulse signal','Impulse response signal');
title('Impulse response signal vs impulse signal');
xlabel('time/s');
ylabel('amp');

%% Plotting b)
NFFT=length(y_unit);

%fft of impulse signal
y_unit_fft=fft(y_unit);
y_unit_fft_m=abs(y_unit_fft);

%fft of impluse response signal
y_response_fft=fft(y_response);
y_response_fft_m=abs(y_response_fft);

%set frequency axis
f_unit=fs*(0:NFFT/2-1)/NFFT;
%% Plotting c)

%show half of the spectrum
y_unit_fft_m=y_unit_fft_m(1:NFFT/2);
y_response_fft_m=y_response_fft_m(1:NFFT/2);

figure;
plot(f_unit,20*log10(y_unit_fft_m));
hold on
plot(f_unit,20*log10(y_response_fft_m));
legend('Impulse signal fft','Impulse response signal fft');
title('Impulse response signal fft vs impulse signal fft');
xlabel('frequency/Hz');
ylabel('amp/dB');
xlim([0 fs/2]);
%ylim([ 0]);
%The fft spectrum of the impulse signal is a constant line while the
%spectrum of the impulse response signal is a faded curve having main lobe and side lobe 
%this is a low pass filter

