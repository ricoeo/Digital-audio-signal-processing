


%% Comb filter a)

fs=44100;   %Hz

%%  Comb filter b)
N=1024;
x_pulse=zeros(1,N);
x_pulse(1)=1;

f_vec= (0:N-1) * fs/N;
comb1=combfilt(x_pulse,-0.5,8);
comb2=combfilt(x_pulse,-0.7,8);
comb3=combfilt(x_pulse,0.5,16);
comb4=combfilt(x_pulse,0.7,16);


figure;
subplot(2,2,1);
%(8,1) (16,-0.5)
stem(0:N-1,comb1);
xlabel('n');
ylabel('h(n)');
xlim([0 128]); 
ylim([-1 1]);
grid on;
title('comb filter with g=-0.5, fold=8 ');
subplot(2,2,2);
%(8,1) (16,-0.7)
stem(0:N-1,comb2);
xlabel('n');
ylabel('h(n)');
xlim([0 128]); 
ylim([-1 1]);
grid on;
title('comb filter with g=-0.7, fold=8 ');
subplot(2,2,3);
%(16,1) (32,-0.5)
stem(0:N-1,comb3);
xlabel('n');
ylabel('h(n)');
xlim([0 128]); 
ylim([-1 1]);
grid on;
title('comb filter with g=0.5, fold=16 ');
subplot(2,2,4);
%(16,1) (32,-0.7)
stem(0:N-1,comb4);
xlabel('n');
ylabel('h(n)');
xlim([0 128]); 
ylim([-1 1]);
grid on;
title('comb filter with g=0.7, fold=16 ');

figure;
plot(f_vec,20*log10(abs(fft(comb1,N))));
hold on;
plot(f_vec,20*log10(abs(fft(comb2,N))));
hold on;
plot(f_vec,20*log10(abs(fft(comb3,N))));
hold on;
plot(f_vec,20*log10(abs(fft(comb4,N))));
legend('g=-0.5, fold=8','g=-0.7, fold=8','g=0.5, fold=16','g=0.7, fold=16')
xlabel('f in Hz');
ylabel('|h(f)|');
xlim([0 fs/2]); ylim([-8 12]); grid on;
title('Magnitude responses of comb filter with different parameters');

%% Allpass filter a)
fs=44100;

%% Allpass filter b)

All1=allfilt(x_pulse,-0.5,8);
All2=allfilt(x_pulse,-0.7,8);
All3=allfilt(x_pulse,-0.5,16);
All4=allfilt(x_pulse,-0.7,16);


figure;
subplot(2,2,1);
stem(0:N-1,All1);
xlabel('n');
ylabel('h(n)');
xlim([0 128]); 
ylim([-1 1]);
grid on;
title('Allpass filter with g=-0.5, fold=8 ');
subplot(2,2,2);
stem(0:N-1,All2);
xlabel('n');
ylabel('h(n)');
xlim([0 128]); 
ylim([-1 1]);
grid on;
title('Allpass filter with g=-0.7, fold=8 ');
subplot(2,2,3);
%(16,1) (32,-0.5)
stem(0:N-1,All3);
xlabel('n');
ylabel('h(n)');
xlim([0 128]); 
ylim([-1 1]);
grid on;
title('Allpass filter with g=-0.5, M=16 ');
subplot(2,2,4);
%(16,1) (32,-0.7)
stem(0:N-1,All4);
xlabel('n');
ylabel('h(n)');
xlim([0 128]); 
ylim([-1 1]);
grid on;
title('Allpass filter with g=-0.7, fold=16 ');

figure;
plot(f_vec,20*log10(abs(fft(All1,N))));
hold on;
plot(f_vec,20*log10(abs(fft(All2,N))));
hold on;
plot(f_vec,20*log10(abs(fft(All3,N))));
hold on;
plot(f_vec,20*log10(abs(fft(All4,N))));
legend('g=-0.5, M=8','g=-0.7, M=8','g=-0.5, M=16','g=-0.7, M=16');
xlabel('f in Hz');
ylabel('|h(f)|');
xlim([0 fs/2]); ylim([-1 1]); grid on;
title('Magnitude responses of allpass filter with different parameters');
%% sound audio 
[audio_vocal,~]=audioread("audio\audio\vocals.wav");
%sound(audio_vocal,Fs);
[audio_funk,Fs]=audioread("audio\audio\funk.wav");
%sound(audio_funk,Fs);
%% Schroeder Reverberator a)
delay=39;   %prime number 30-45
g=-0.85;

audio_filter_comb1=combfilt(audio_vocal,-0.45,delay);
audio_filter_comb2=combfilt(audio_filter_comb1,g,20);
audio_filter_comb3=combfilt(audio_filter_comb2,-0.3,15);
audio_filter_comb4=combfilt(audio_filter_comb3,0.3,10);

sound(audio_filter_comb4,Fs);


%% Schroeder Reverberator a)
delay=5;
g=-0.7;

audio_filter_all1=allfilt(audio_vocal,g,delay);
audio_filter_all2=allfilt(audio_filter_all1,g,3);

sound(audio_filter_all2,Fs);

%% Schroeder Reverberator c)
%Building filter structure
my_Schroeder_unit_response=my_Schroeder_reverberator(x_pulse);
figure;
stem(0:N-1,my_Schroeder_unit_response);
xlabel('n');
ylabel('h(n)');
xlim([0 128]); 
ylim([-1 1]);
grid on;
title('Impulse responese of my schroeder reverberator');
%% 
figure;
my_Schroeder_unit_fft=abs(fft(my_Schroeder_unit_response,N));
semilogx(f_vec,my_Schroeder_unit_fft);
grid on;
xlim([0 fs/2]); ylim([-2 14]);
title("Frequency response of my schroeder reverberator");

%% 

sound(my_Schroeder_reverberator(audio_vocal),Fs);

%% 

sound(my_Schroeder_reverberator(audio_funk),Fs);

%% 

NFFT=length(audio_vocal);
f_unit= (0:floor(NFFT/2)-1) * fs/NFFT;
audio_fft=abs(fft(audio_vocal,NFFT));
audio_fft=audio_fft(1:floor(NFFT/2));
figure;
semilogx(f_unit,20*log10(audio_fft));
title('Frequency response of audio')
xlabel('f'); ylabel('|Y(f)| / dB');

%% functions

function [y]=combfilt(x,g,M)
    % Comb filter
    b = zeros(1,M+1);    	% Numerator coefficients
    b(length(b))= 1;
    a = zeros(1,M+1);     % Denominator coefficients
    a(1)= 1;
    a(length(a))= -g;
    y= filter(b,a,x);    % Impulse response
end

function [y]=allfilt(x,g,M)
    % Comb filter
    b = zeros(1,M+1);    	% Numerator coefficients
    b(1)=-g;
    b(length(b))= 1;
    a = zeros(1,M+1);     % Denominator coefficients
    a(1)= 1;
    a(length(a))= -g;
    y= filter(b,a,x);    % Impulse response
end

function [y]=my_Schroeder_reverberator(x)
    audio_filter1=combfilt(x,-0.45,39);
    audio_filter2=combfilt(audio_filter1,-0.85,20);
    audio_filter3=combfilt(audio_filter2,-0.3,15);
    audio_filter4=combfilt(audio_filter3,0.3,10);
    audio_filter5=allfilt(audio_filter4,-0.7,13);
    audio_filter6=allfilt(audio_filter5,0.3,7);
    audio_filter7=allfilt(audio_filter6,-0.3,5);
    y=allfilt(audio_filter7,0.6,3);
end
