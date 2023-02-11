%% Shelving Filter a)
%high frequency shelving filter with cut

fs=44100;   %unit Hz
fc=3000;    %unit Hz
G= 15;      %unit dB

T=1/fs;
N=4096;
n= [1,zeros(1,N)]; 	% dirac impulse
f_vec= (0:N-1) * fs/N;

V0b= 10^(-G/20);
H0b= V0b - 1;
a_c= (1-V0b*tan(2*pi*fc*T/2)) /(V0b* tan(2*pi*fc*T/2)+1);

% All-pass
Bc = [-a_c, 1];
Ac = [1, -a_c];
y_all = filter(Bc,Ac,n);
y_c= n + H0b/2 * (n-y_all);
%% Shelving filter b)
figure;
plot(f_vec,20*log10(abs(fft(y_c,N))));
title('Magnitude response of high-freuency shelving (cut) output')
xlim([0 fs/2]); 
ylim([-G G]); 
grid on;
xlabel('f in Hz'); ylabel('|Y_C(f)| in dB');

%% Shelving filter c)
figure;
semilogx(f_vec,20*log10(abs(fft(y_c,N))));
title('Magnitude response of high-freuency shelving (cut) output')
xlim([0 fs/2]); ylim([-G G]); grid on;
xlabel('f in Hz'); ylabel('|Y_C(f)| in dB');

%% Shelving filter d)
test=highshelvingfilt(n,-15,3000,44100);
figure;
semilogx(f_vec,20*log10(abs(fft(test,N))));
%HF shelving filter with cut case


%% Shelving filter e)
figure;
subplot(1,2,1);
test1=highshelvingfilt(n,-12,3000,44100);
semilogx(f_vec,20*log10(abs(fft(test1,N))));
hold on
test2=highshelvingfilt(n,-15,3000,44100);
semilogx(f_vec,20*log10(abs(fft(test2,N))));
hold on
test3=highshelvingfilt(n,-18,3000,44100);
semilogx(f_vec,20*log10(abs(fft(test3,N))));
hold off;
legend('G=-12dB','G=-15dB','G=-18dB');
title("HF shelving filters with different gain");
xlabel('frequency');
ylabel('|Y_C(f)| in dB');
xlim([0 fs/2]);


subplot(1,2,2);
test4=highshelvingfilt(n,-12,1000,44100);
semilogx(f_vec,20*log10(abs(fft(test4,N))));
hold on
test5=highshelvingfilt(n,-12,2000,44100);
semilogx(f_vec,20*log10(abs(fft(test5,N))));
hold on
test6=highshelvingfilt(n,-12,3000,44100);
semilogx(f_vec,20*log10(abs(fft(test6,N))));
hold off;
legend('fc=1kHz','fc=2kHz','fc=3kHz');
title("HF shelving filters with different cut-off frequency");
xlabel('frequency');
ylabel('|Y_C(f)| in dB');
xlim([0 fs/2]);
%% Peak filter a)

fs=44100;
G=6;
Q=4;
fc=1000;

T=1/fs;
fb=fc/Q;
V0p=10^(G/20);
H0p=V0p-1;
b=cos(2*pi*fc*T);
a_p=(1-tan(2*pi*fb*T/2))/(1+tan(2*pi*fb*T/2));
%all-pass
B_p=[a_p,-b*(1+a_p),1];
A_p=[1,-b*(1+a_p),a_p];
y_p_all=filter(B_p,A_p,n);

y_p= n + H0p/2 * (n-y_p_all);
%% Peak filter b)c)
N=4096;
figure;
semilogx(f_vec,20*log10(abs(fft(y_p,N))));
xlim([0 fs/2]); 
ylim([0 G]); 
grid on;
xlabel('f in Hz'); ylabel('|Y_C(f)| in dB');
title('Magnitude response of peak filter');
%% Peak filter d)
test=peakfilt(n,6,1000,4,44100);
figure;
semilogx(f_vec,20*log10(abs(fft(test,N))));

%% Peak filter e)
figure;
subplot(1,2,1);
test1=peakfilt(n,6,1000,2,44100);
semilogx(f_vec,20*log10(abs(fft(test1,N))));
hold on
test2=peakfilt(n,6,1000,3,44100);
semilogx(f_vec,20*log10(abs(fft(test2,N))));
hold on
test3=peakfilt(n,6,1000,4,44100);
semilogx(f_vec,20*log10(abs(fft(test3,N))));
hold off;
legend('Q=2','Q=3','Q=4');
xlim([0 fs/2]); 
grid on;
xlabel('f in Hz'); ylabel('|Y_C(f)| in dB');
title('peak filter with different Q');

subplot(1,2,2);
test4=peakfilt(n,-6,1000,4,44100);
semilogx(f_vec,20*log10(abs(fft(test4,N))));
hold on
test5=peakfilt(n,-6,1200,4,44100);
semilogx(f_vec,20*log10(abs(fft(test5,N))));
hold on
test6=peakfilt(n,-6,800,4,44100);
semilogx(f_vec,20*log10(abs(fft(test6,N))));
hold off;
legend('fc=1000Hz','fc=1200Hz','fc=800Hz');
xlim([0 fs/2]); 
grid on;
xlabel('f in Hz'); ylabel('|Y_C(f)| in dB');
title('peak filter with different cut-off frequency');

%% Parametric Equaliser a)
[audio,Fs]=audioread("audio\audio\drumloop.wav");
f_unit= (0:176400-1) * fs/176400;
figure;
semilogx(f_unit,20*log10(abs(fft(audio,176400))));
hold on
y_temp1=highshelvingfilt(audio,-2,400,Fs);
y_temp2=peakfilt(y_temp1,6,1000,2,Fs);
y_temp3=peakfilt(y_temp2,6,1000,2,Fs);
y=highshelvingfilt(y_temp3,-12,5000,Fs); %dampen
sound(y,fs);
semilogx(f_unit,20*log10(abs(fft(y,176400))));
hold off;
legend('audio','filtered audio');
xlabel('f');
ylabel('|Y_C(f)| in dB');
%% Parametric Equaliser b)
[audio_pop,Fs]=audioread("audio\audio\popsong.wav");
N_pop=length(audio_pop);
f_unit_pop=(0:N_pop-1)*fs/N_pop;
n_pop=0:N_pop-1;
n_pop=n_pop';
f_noise=440;
%create noise signal
x_noise=0.5*sin(n_pop*2*pi*f_noise/fs);
audio_pop_rev=audio_pop+x_noise;
%sound(audio_pop_rev,fs)
Q_select=6;
gain_select=30;
audio_pop_fil=peakfilt(audio_pop_rev,-gain_select,f_noise,Q_select,fs);
figure;
semilogx(f_unit_pop,20*log10(abs(fft(audio_pop_rev,N_pop))));
hold on;
semilogx(f_unit_pop,20*log10(abs(fft(audio_pop_fil,N_pop))));
hold off;
legend('audio','filter');
%after the filter, the noise is nearly filtered out from the original signal
sound(audio_pop_fil,fs);

%% Parametric Equaliser c)

%Drum low frequency enhance
Drum_filt=myEQ1(audio);
sound(Drum_filt,fs);
%% 

%Pop song more warm
Pop_filt=myEQ2(audio_pop);
sound(Pop_filt,fs);
%% function
function [y] = myEQ1(x)
    temp1=peakfilt(x,6,100,3,44100);
    y=highshelvingfilt(temp1,-6,3000,44100);
end
function [y] = myEQ2(x)
    temp1=lowshelvingfilt(x,-10,250,44100);
    temp2=peakfilt(temp1,-6,500,2,44100);
    temp3=peakfilt(temp2,3,240,7,44100);
    y=highshelvingfilt(temp3,-10,4000,44100);
end
function [y] =highshelvingfilt(x, gain, fc,fs)
    V0= 10^(gain/20);  
    H0= V0 - 1;
    a=(V0*tan(pi*fc/fs)-1)/(1+V0*tan(pi*fc/fs));
    y_all_temp=filter([a,1],[1,a],x);
    y= x + H0/2 * (x-y_all_temp);
end
function [y] =lowshelvingfilt(x, gain, fc,fs)
    V0= 10^(gain/20);   
    H0= V0 - 1;
    a=(V0*tan(pi*fc/fs)-1)/(1+V0*tan(pi*fc/fs));
    y_all_temp=filter([a,1],[1,a],x);
    y= x + H0/2 * (x+y_all_temp);
end
function [y] = peakfilt(x, gain, fc, Q,fs)
    fb=fc/Q;
    V0=10^(gain/20);
    H0=V0-1;
    b=cos(2*pi*fc/fs);
    a=(1-tan(pi*fb/fs))/(1+tan(pi*fb/fs));
    %all-pass
    B=[a,-b*(1+a),1];
    A=[1,-b*(1+a),a];
    y_all=filter(B,A,x);
    y= x + H0/2 * (x-y_all);
end