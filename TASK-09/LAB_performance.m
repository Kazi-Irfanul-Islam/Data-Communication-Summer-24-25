% ======================================================
% Lab Performance Task - Frequency Division Multiplexing (4 Signals)
% ======================================================

clc; clear; close all;

%% User Defined Parameters
F = 2; % Example, replace with actual value
G = 3; % Example, replace with actual value

fs = 4001; % Sampling Frequency
t = 0:1/fs:1-1/fs; % Time axis

%% Message Signal Generation
am1 = F + 2; fm1 = G + 1;
am2 = F + 5; fm2 = G + 2;
am3 = F + 8; fm3 = G + 3;
am4 = F + 11; fm4 = G + 4;

mt1 = am1*cos(2*pi*fm1*t);
mt2 = am2*cos(2*pi*fm2*t);
mt3 = am3*cos(2*pi*fm3*t);
mt4 = am4*cos(2*pi*fm4*t);

%% Carrier Signal Generation
% Choose carriers within 50-250 Hz range and spaced apart
fc1 = 60; fc2 = 110; fc3 = 170; fc4 = 230; 
Cm = 1; % Carrier amplitude

c1 = Cm*cos(2*pi*fc1*t);
c2 = Cm*cos(2*pi*fc2*t);
c3 = Cm*cos(2*pi*fc3*t);
c4 = Cm*cos(2*pi*fc4*t);

%% Amplitude Modulation (DSB-SC)
s1 = mt1 .* c1;
s2 = mt2 .* c2;
s3 = mt3 .* c3;
s4 = mt4 .* c4;

%% Composite Signal (FDM)
x = s1 + s2 + s3 + s4;

%% Plot Time Domain of Message Signals
figure
subplot(4,1,1); plot(t, mt1); title('Message Signal 1'); xlabel('Time'); ylabel('Amplitude'); ylim([-am1 am1])
subplot(4,1,2); plot(t, mt2); title('Message Signal 2'); xlabel('Time'); ylabel('Amplitude'); ylim([-am2 am2])
subplot(4,1,3); plot(t, mt3); title('Message Signal 3'); xlabel('Time'); ylabel('Amplitude'); ylim([-am3 am3])
subplot(4,1,4); plot(t, mt4); title('Message Signal 4'); xlabel('Time'); ylabel('Amplitude'); ylim([-am4 am4])

%% Frequency Axis
f = fs/2*linspace(-1,1,fs);

%% FFT of Composite Signal
X = abs(fftshift(fft(x)))/(fs/2);
figure
stem(f,X); title('Composite Signal Spectrum'); xlabel('Frequency'); ylabel('Amplitude'); xlim([-270 270])

%% Bandpass Filtering (Demultiplexing)
[num1, den1] = butter(5, [(fc1-fm1-5)/(fs/2),(fc1+fm1+5)/(fs/2)]);
[num2, den2] = butter(5, [(fc2-fm2-5)/(fs/2),(fc2+fm2+5)/(fs/2)]);
[num3, den3] = butter(5, [(fc3-fm3-5)/(fs/2),(fc3+fm3+5)/(fs/2)]);
[num4, den4] = butter(5, [(fc4-fm4-5)/(fs/2),(fc4+fm4+5)/(fs/2)]);

bpf1 = filter(num1, den1, x);
bpf2 = filter(num2, den2, x);
bpf3 = filter(num3, den3, x);
bpf4 = filter(num4, den4, x);

%% Mixing with Carrier for Demodulation
z1 = 2 * bpf1 .* c1;
z2 = 2 * bpf2 .* c2;
z3 = 2 * bpf3 .* c3;
z4 = 2 * bpf4 .* c4;

%% Lowpass Filtering to Recover Messages
[num_lp1, den_lp1] = butter(5, (fm1+3)/(fs/2));
[num_lp2, den_lp2] = butter(5, (fm2+3)/(fs/2));
[num_lp3, den_lp3] = butter(5, (fm3+3)/(fs/2));
[num_lp4, den_lp4] = butter(5, (fm4+3)/(fs/2));

rec1 = filter(num_lp1, den_lp1, z1);
rec2 = filter(num_lp2, den_lp2, z2);
rec3 = filter(num_lp3, den_lp3, z3);
rec4 = filter(num_lp4, den_lp4, z4);

%% Plot Received Signals (Time Domain)
figure
subplot(4,1,1); plot(t, rec1); title('Recovered Message 1'); xlabel('Time'); ylabel('Amplitude'); ylim([-am1 am1])
subplot(4,1,2); plot(t, rec2); title('Recovered Message 2'); xlabel('Time'); ylabel('Amplitude'); ylim([-am2 am2])
subplot(4,1,3); plot(t, rec3); title('Recovered Message 3'); xlabel('Time'); ylabel('Amplitude'); ylim([-am3 am3])
subplot(4,1,4); plot(t, rec4); title('Recovered Message 4'); xlabel('Time'); ylabel('Amplitude'); ylim([-am4 am4])

%% Plot Frequency Domain of Received Signals
R1 = abs(fftshift(fft(rec1)))/(fs/2);
R2 = abs(fftshift(fft(rec2)))/(fs/2);
R3 = abs(fftshift(fft(rec3)))/(fs/2);
R4 = abs(fftshift(fft(rec4)))/(fs/2);

figure
subplot(4,1,1); stem(f, R1); title('Recovered Message 1 Spectrum'); xlabel('Frequency'); ylabel('Amplitude'); xlim([-10 10])
subplot(4,1,2); stem(f, R2); title('Recovered Message 2 Spectrum'); xlabel('Frequency'); ylabel('Amplitude'); xlim([-10 10])
subplot(4,1,3); stem(f, R3); title('Recovered Message 3 Spectrum'); xlabel('Frequency'); ylabel('Amplitude'); xlim([-10 10])
subplot(4,1,4); stem(f, R4); title('Recovered Message 4 Spectrum'); xlabel('Frequency'); ylabel('Amplitude'); xlim([-10 10])
