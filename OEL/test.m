

clear;
clc; 
close all;

%% Task 1: Generate Binary Vectors
surname1 = 'IRFAN';
surname2 = 'SIUM';

bitsChar1 = dec2bin(uint8(surname1), 8).'; 
x1 = bitsChar1(:).' - '0';                   

bitsChar2 = dec2bin(uint8(surname2), 8).';
x2 = bitsChar2(:).' - '0';

disp('Binary vector x1 (IRFAN):'); disp(x1);
disp('Binary vector x2 (SIUM):'); disp(x2);

Nb1 = numel(x1);      
Nb2 = numel(x2);      


%% Task 2: NRZ Waveform
Tb  = 1e-3;       
Fs  = 100e3;       
Ns  = round(Fs*Tb);
A_nrz = 1;         

nrz1 = A_nrz * kron(x1, ones(1, Ns));
nrz2 = A_nrz * kron(x2, ones(1, Ns));

t1 = (0:numel(nrz1)-1)/Fs;
t2 = (0:numel(nrz2)-1)/Fs;

% Plot NRZ
figure(1); clf;
subplot(2,2,1); plot(t1*1e3, nrz1,'LineWidth',1.2); grid on; ylim([-0.2 1.2]);
xlabel('Time (ms)'); ylabel('Amplitude'); title('x1 (IRFAN) - NRZ');

subplot(2,2,2); plot(t2*1e3, nrz2,'LineWidth',1.2); grid on; ylim([-0.2 1.2]);
xlabel('Time (ms)'); ylabel('Amplitude'); title('x2 (SIUM) - NRZ');

%% Task 3: BASK Modulation
fc1 = 2000;      % Carrier 1 = 2 kHz
fc2 = 4000;      % Carrier 2 = 4 kHz
Ac  = 1;         

c1 = cos(2*pi*fc1*t1);
c2 = cos(2*pi*fc2*t2);

s1_fc2kHz = Ac * nrz1 .* c1;
s2_fc4kHz = Ac * nrz2 .* c2;

subplot(2,2,3); plot(t1, s1_fc2kHz,'LineWidth',1.1); grid on;
xlabel('Time (s)'); ylabel('Amplitude'); title('BASK x1 @ 2 kHz');

subplot(2,2,4); plot(t2, s2_fc4kHz,'LineWidth',1.1); grid on;
xlabel('Time (s)'); ylabel('Amplitude'); title('BASK x2 @ 4 kHz');

%% Task 4: Composite Signal + Gausian Noise
len1 = length(s1_fc2kHz); len2 = length(s2_fc4kHz); L = max(len1,len2);

if len1<L, s1_fc2kHz = [s1_fc2kHz zeros(1,L-len1)]; end
if len2<L, s2_fc4kHz = [s2_fc4kHz zeros(1,L-len2)]; end

composite_clean = s1_fc2kHz + s2_fc4kHz;

SNR_dB = 15;
P_sig = mean(composite_clean.^2);
sigma2 = P_sig / (10^(SNR_dB/10));
noise = sqrt(sigma2) * randn(size(composite_clean));
composite_noisy = composite_clean + noise;

figure(2); clf;
subplot(3,1,1); plot((0:L-1)/Fs, composite_noisy,'LineWidth',1); grid on;
xlabel('Time (s)'); ylabel('Amplitude'); title('Composite Signal with AWGN');

%% Task 5: Demultiplexing & Demodulation 
Ns = round(Fs*Tb);
lp_len = round(0.75*Ns); if mod(lp_len,2)==0, lp_len=lp_len+1; end
lpf = ones(1,lp_len)/lp_len;

composite = composite_noisy;
tp = (0:length(composite)-1)/Fs;

% Demodulate x1
bb1 = composite .* cos(2*pi*fc1*tp);
env1 = conv(bb1,lpf,'same');
delay = floor((lp_len-1)/2);
env1 = [env1(delay+1:end) zeros(1,delay)];

% Demodulate x2
bb2 = composite .* cos(2*pi*fc2*tp);
env2 = conv(bb2,lpf,'same');
env2 = [env2(delay+1:end) zeros(1,delay)];

% Sample bit centers
centers = round(Ns/2 : Ns : length(env1)-Ns/2);
samp1 = env1(centers(1:Nb1)); samp2 = env2(centers(1:Nb2));
x1_hat_raw = samp1 > 0.5*max(samp1);
x2_hat_raw = samp2 > 0.5*max(samp2);

% --- Alignmentfor x1 ---
best_shift1 = 0; best_corr1 = -Inf;
for shift=0:7
    trial = x1_hat_raw(1+shift:end); trial = trial(1:min(Nb1,numel(trial)));
    r = sum(trial == x1(1:length(trial)));
    if r>best_corr1, best_corr1=r; best_shift1=shift; end
end
x1_hat = x1_hat_raw(1+best_shift1 : min(Nb1+best_shift1, numel(x1_hat_raw)));

% --- Alignment for x2 ---
best_shift2 = 0; best_corr2 = -Inf;
for shift=0:7
    trial = x2_hat_raw(1+shift:end); trial = trial(1:min(Nb2,numel(trial)));
    r = sum(trial == x2(1:length(trial)));
    if r>best_corr2, best_corr2=r; best_shift2=shift; end
end
x2_hat = x2_hat_raw(1+best_shift2 : min(Nb2+best_shift2, numel(x2_hat_raw)));

% Convert to NRZ waveform
nrz1_rec = A_nrz * kron(x1_hat, ones(1,Ns));
nrz2_rec = A_nrz * kron(x2_hat, ones(1,Ns));
t_rec1 = (0:length(nrz1_rec)-1)/Fs;
t_rec2 = (0:length(nrz2_rec)-1)/Fs;

subplot(3,1,2);
stairs(t_rec1*1e3, nrz1_rec,'LineWidth',1.2); grid on; ylim([-0.2 1.2]);
xlabel('Time (ms)'); ylabel('Amplitude'); title('Recovered NRZ x1 (IRFAN)');

subplot(3,1,3);
stairs(t_rec2*1e3, nrz2_rec,'LineWidth',1.2); grid on; ylim([-0.2 1.2]);
xlabel('Time (ms)'); ylabel('Amplitude'); title('Recovered NRZ x2 (SIUM)');

%% Task 6: ASCII Decoding & Verification
bits2chars = @(b) char(bi2de(reshape(b(1:8*floor(numel(b)/8)),8,[]).','left-msb')).';

msg1     = bits2chars(x1);
msg2     = bits2chars(x2);
msg1_hat = bits2chars(x1_hat);
msg2_hat = bits2chars(x2_hat);

fprintf('\nOriginal x1 (IRFAN): %s\n', msg1);
fprintf('Recovered x1: %s\n', msg1_hat);

fprintf('\nOriginal x2 (SIUM): %s\n', msg2);
fprintf('Recovered x2: %s\n', msg2_hat);
