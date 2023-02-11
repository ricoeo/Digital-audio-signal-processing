%% Quantisation a)
w=16;
Q=2^(-(w-1));       %quant step

%% Quantisation b)
n=-1:1/1000:1;
x=3*Q*n;    %x(n)
xQ1=floor(x/Q+0.5)*Q;
xQ2=floor(x/Q)*Q;
%characteristic curves
figure;
plot(n,x);
hold on
plot(n,xQ1);
hold on
plot(n,xQ2);

legend('Input Data','Quantisation with rounding','Quantisation with truncation');
title('characteristic curves for two Quantization types');
xlabel('n');
ylabel('amp');
ylim([-3*Q 3*Q]);

%error signals
e1=xQ1-x;
e2=xQ2-x;
figure;
plot(n,e1);
hold on
plot(n,e2);
legend('Error with rounding','Error with truncation');
title('Error curves for two Quantization types');
xlabel('n');
ylabel('amp');

%% Quantisation c)

N=1024;
fs=44100;   %unit Hz
n=0:N-1;

x = Q* sin(2*pi* 64/N *n);  %input x(n)

XQ3_qt=floor(x/Q)*Q;
XQ3_qr=floor(x/Q+0.5)*Q;
n_show=0:99;

figure;
plot(XQ3_qt(n_show+1));
hold on;
plot(XQ3_qr(n_show+1));
hold on;
plot(x(n_show+1));
hold off;
legend('quantization with truncation','quantization with rounding','sin signal');
title('Quantisation of input signal x(n)');
xlabel('n'); ylabel('x_{Q}(n)');
%% 

NFFT=N;
x_withditter_fft=fft(XQ3_qt,NFFT);
XQ3_qr_fft=fft(XQ3_qr,NFFT);

x_withditter_m=abs(x_withditter_fft);
XQ3_qr_fft_m=abs(XQ3_qr_fft);

x_withditter_m=x_withditter_m(1:NFFT/2);
XQ3_qr_fft_m=XQ3_qr_fft_m(1:NFFT/2);
f_unit=fs*(0:NFFT/2-1)/NFFT;
figure;
subplot(1,2,1);
plot(f_unit,20*log10(x_withditter_m));
xlabel('f'); ylabel('|X_{Q}(f)| / dB');
subplot(1,2,2);
plot(f_unit,20*log10(abs(XQ3_qr_fft_m)+eps));    %nothing?
xlabel('f/Hz'); ylabel('|X_{Q}(f)| / dB');

%% Quantisation d)
e_qt=x-XQ3_qt;
e_qr=x-XQ3_qr;
p_Et= 1/length(e_qt) * sum(e_qt.^2);    % Error power
p_Er= 1/length(e_qr) * sum(e_qr.^2);
p_X= 1/length(x) * sum(x.^2);     % Signal power
SNRt= 10* log10(p_X/p_Et);
SNRr= 10* log10(p_X/p_Er);

fprintf('the quantization error power for Truncation method is %d w. The SNR is %f dB. \n',p_Et,SNRt);
fprintf('the quantization error power for Rounding method is %d w. The SNR is %f dB.\n',p_Er,SNRr);

%% Quantisation e)
a = [1,0.5, 0.1, Q];
x_re1 = a(1)* sin(2*pi* 64/N *n);
x_re2 = a(2)* sin(2*pi* 64/N *n);
x_re3 = a(3)* sin(2*pi* 64/N *n);
%x_re4 = a(4)* sin(2*pi* 64/N *n);

x_re1_qr=floor(x_re1/Q+0.5)*Q;
x_re2_qr=floor(x_re2/Q+0.5)*Q;
x_re3_qr=floor(x_re3/Q+0.5)*Q;
x_re1_qt=floor(x_re1/Q)*Q;
x_re2_qt=floor(x_re2/Q)*Q;
x_re3_qt=floor(x_re3/Q)*Q;

e_re1_qr=x_re1_qr-x_re1;
e_re2_qr=x_re2_qr-x_re2;
e_re3_qr=x_re3_qr-x_re3;
e_re1_qt=x_re1_qt-x_re1;
e_re2_qt=x_re2_qt-x_re2;
e_re3_qt=x_re3_qt-x_re3;

p_Et1= 1/length(e_re1_qt) * sum(e_re1_qt.^2);
p_Er1= 1/length(e_re1_qr) * sum(e_re1_qr.^2);
p_X1= 1/length(x_re1) * sum(x_re1);    
SNRt1= 10* log10(p_X1/p_Et1);
SNRr1= 10* log10(p_X1/p_Er1);

p_Et2= 1/length(e_re2_qt) * sum(e_re2_qt.^2);
p_Er2= 1/length(e_re2_qr) * sum(e_re2_qr.^2);
p_X2= 1/length(x_re2) * sum(x_re2);    
SNRt2= 10* log10(p_X2/p_Et2);
SNRr2= 10* log10(p_X2/p_Er2);

p_Et3= 1/length(e_re3_qt) * sum(e_re3_qt.^2);
p_Er3= 1/length(e_re3_qr) * sum(e_re3_qr.^2);
p_X3= 1/length(x_re3) * sum(x_re3);    
SNRt3= 10* log10(p_X3/p_Et3);
SNRr3= 10* log10(p_X3/p_Er3);

p_ET=[p_Et1;p_Et2;p_Et3;p_Et];
p_ER=[p_Er1;p_Er2;p_Er3;p_Er];
SNRT=[SNRt1;SNRt2;SNRt3;SNRt];
SNRR=[SNRr1;SNRr2;SNRr3;SNRr];

T1=table(p_ET,SNRT,'VariableNames',{'Truncation method:quantization error','SNR'},'RowNames',{'a=1';'a=0.5';'a=0.1';'a=Q'});
T2=table(p_ER,SNRR,'VariableNames',{'Rounding method:quantization error','SNR'},'RowNames',{'a=1';'a=0.5';'a=0.1';'a=Q'});

%% Quantisation f)
[audio,Fs]=audioread("audio\audio\drumloop.wav");
W=[1,2,3,4];
[audio_re1,p1,SNR1]=sig_process(audio,W(1),"round");
[audio_re2,p2,SNR2]=sig_process(audio,W(2),"round");
[audio_re3,p3,SNR3]=sig_process(audio,W(3),"round");
[audio_re4,p4,SNR4]=sig_process(audio,W(4),"round");

[audio_re5,p5,SNR5]=sig_process(audio,W(1),"trun");
[audio_re6,p6,SNR6]=sig_process(audio,W(2),"trun");
[audio_re7,p7,SNR7]=sig_process(audio,W(3),"trun");
[audio_re8,p8,SNR8]=sig_process(audio,W(4),"trun");
%sound(audio_re1,Fs);
%sound(audio_re3,Fs);
%sound(audio_re5,Fs);
%sound(audio_re4,Fs);
%sound(audio_re8,Fs);

[audio_re9,p9,SNR9]=sig_process(audio,128,"round");
sound(audio_re9,Fs);


%% Dithering a)
NFFT        = 1024;                      

% Rectangular dither
r1          = 2 * rand(1,NFFT) - 1;
d_rect      = r1 * Q/2; 

%amplitude distribution
f_vec               = (0:(NFFT-1)) / NFFT * fs;
[n_hist1, d_hist1]  = hist(d_rect, 11);
%frequency spectrum
D_rect              = fft(d_rect,NFFT);

figure;
subplot(1,2,1);
bar(d_hist1/Q, 100*n_hist1/sum(n_hist1));
xlim([-1 1]); ylim([0 12]);
xlabel('d_{rect}/Q'); ylabel('Relative frequency in %');
title('amplitude distribution');
subplot(1,2,2);
plot(f_vec, 20*log10(abs(D_rect)));
xlim([0 f_s/2]); ylim([-110 -60]); 
title('frequency spectrum');
xlabel('f'); 
ylabel('|D_{rect}(f)| /dB');

%% Dithering b)
x_withditter=x+d_rect;

x_withditter_qr=floor(x_withditter/Q+0.5)*Q;
figure;
plot(x_withditter_qr(n_show+1));
title('Quantisation of x(n) with ditter');
xlabel('n'); ylabel('x_{Q}(n)');
%% 
x_withditter_fft=fft(x_withditter,NFFT);
x_withditter_m=abs(x_withditter_fft);
x_withditter_m=x_withditter_m(1:NFFT/2);

figure;
plot(f_unit,20*log10(x_withditter_m));
title('Frequency response of quantized x(n)')
xlabel('f'); ylabel('|X_{Q}(f)| / dB');

%% Dithering c)
e_withditter=x_withditter_qr-x;
p_ED=1/length(e_withditter) * sum(e_withditter.^2);
p_XD=1/length(x_withditter) * sum(x_withditter.^2);
SNR_ED=10* log10(p_XD/p_ED);
fprintf('The error power without dittering is %d, with dittering is %d \n',p_Er,p_ED);
fprintf('The SNR without dittering is %f dB, with dittering is %f dB \n',SNRr,SNR_ED);

%% Dithering d)
[audio_pop,Fs]=audioread("audio\audio\popsong.wav");
W=[1,2,3,4];
sound(audio_pop,Fs);
%% 
%audio1=sig_process_ditter(audio_pop,W(1),"on");
%audio1=sig_process_ditter(audio_pop,W(4),"on");
%audio1=sig_process_ditter(audio_pop,64,"on");
audio1=sig_process_ditter(audio_pop,W(4),"off");
sound(audio1,Fs);
%% function sig_process
function [audio_re,p,SNR]=sig_process(original_audio,w_opt,method)
    q_opt=2^(-(w_opt-1));
    if(method=="trun")
        audio_re=floor(original_audio/q_opt)*q_opt;
    elseif(method=="round")
        audio_re=floor(original_audio/q_opt+0.5)*q_opt;
    end
    p=1/length((original_audio-audio_re)) * sum((original_audio-audio_re).^2);
    p_sig=1/length(original_audio) * sum(original_audio);
    SNR=10* log10(p_sig/p);
end

function [audio_re]=sig_process_ditter(original_audio,w_opt,switch_ditter)
    q_opt=2^(-(w_opt-1));
    ditter=(2 * rand(length(original_audio),1) - 1)*q_opt/2;
    if(switch_ditter=="on")
        audio_re=floor((original_audio+ditter)/q_opt+0.5)*q_opt;
    elseif(switch_ditter=="off")
        audio_re=floor(original_audio/q_opt+0.5)*q_opt;
    end
end
