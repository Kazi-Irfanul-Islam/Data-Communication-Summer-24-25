% ======================================================
% Lab Performance Task - QASK
% Transmit and Receive 'Data Comm is fun'
% ======================================================

clc; clear; close all;

%% (a) Transmit text message as analog signal using QASK
msg = 'Data Comm is fun';
disp('Original Message:');
disp(msg);

% Convert to ASCII binary
binData = dec2bin(msg,8)'; 
binData = binData(:)'; % row vector
binData = binData - '0'; % convert '0'/'1' to numeric

disp('Binary Data:');
disp(binData);

% Parameters
M = 16;                    % QASK (16-QAM)
k = log2(M);               % Bits per symbol
binData = [binData zeros(1,mod(-length(binData),k))]; % pad to multiple of k

% Digital Signal Representation (NRZ)
digitalSig = [];
for n=1:length(binData)
    if binData(n)==1
        digitalSig = [digitalSig ones(1,100)];
    else
        digitalSig = [digitalSig zeros(1,100)];
    end
end

% Plot binary digital signal
figure;
subplot(3,1,1);
plot(digitalSig,'LineWidth',1.5);
ylim([-0.5 1.5]); grid on;
title('Digital Signal (Binary Data Representation)');
xlabel('Samples'); ylabel('Amplitude');

% QASK Modulation
dataSymbols = bi2de(reshape(binData,k,length(binData)/k).','left-msb');
modSig = qammod(dataSymbols,M,'UnitAveragePower',true);

% ===== Upsample to create continuous analog waveform =====
upsampleFactor = 50; % samples per symbol
modSigUp = upsample(modSig,upsampleFactor);
filt = rcosdesign(0.35,6,upsampleFactor); % pulse shaping filter
txWaveform = conv(modSigUp,filt,'same');

% Plot I/Q constellation and analog waveform
subplot(3,1,2);
plot(real(modSig), imag(modSig),'o');
title('QASK Constellation (16-QAM)');
xlabel('In-phase'); ylabel('Quadrature'); grid on;

subplot(3,1,3);
N = min(200,length(txWaveform));   % take available samples
plot(real(txWaveform(1:N)),'LineWidth',1.5);
title('Analog QASK Waveform (Time Domain)');
xlabel('Samples'); ylabel('Amplitude'); grid on;

%% (b) Channel with SNR = 30 dB
rxSig = awgn(txWaveform,30,'measured');

figure;
plot(real(rxSig(1:200)),'LineWidth',1.5);
title('Received QASK Signal with Noise (SNR=30 dB)');
xlabel('Samples'); ylabel('Amplitude'); grid on;

%% (c) Demodulation and Text Recovery
% Match filter (same rcosdesign)
rxFilt = conv(rxSig,filt,'same');

% Downsample back to symbols
rxDown = downsample(rxFilt,upsampleFactor);

% QASK Demodulation
rxDataSymbols = qamdemod(rxDown,M,'UnitAveragePower',true);
rxBits = de2bi(rxDataSymbols,k,'left-msb')';
rxBits = rxBits(:)';

% Convert bits back to characters
rxBits = rxBits(1:length(msg)*8); % remove padding
rxChars = char(bin2dec(reshape(char(rxBits+'0'),8,[]).')).';

disp('Recovered Text:');
disp(rxChars);

% Show only "Data Comm"
disp('Recovered Substring:');
disp(rxChars(1:9));
