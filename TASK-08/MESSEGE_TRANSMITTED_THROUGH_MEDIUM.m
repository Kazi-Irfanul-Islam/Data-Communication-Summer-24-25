disp('**********')
disp(' Message transmitted through a Transmission medium');
disp('**********')

% Channel Noise %
t4 = bp/99:bp/99:bp*length(x);
Rec = awgn(m,10); % Additive White Gaussian Noise
subplot(4,1,3);
plot(t4,Rec);
axis([0 bp*length(x) -6 6]);
xlabel('time(sec)');
ylabel('amplitude(volt)');
title('Received signal at Receiver');

% ========== Binary ASK Demodulation ==========
mn = [];
for n = ss:ss:length(Rec)
    t = bp/99:bp/99:bp;
    y = cos(2*pi*f*t);        % carrier signal
    mm = y .* Rec((n-(ss-1)):n); 
    t5 = bp/99:bp/99:bp;
    z = trapz(t5,mm);         % integration 
    zz = round((2*z/bp)); 
    
    if (zz > 2.5)             % logic level threshold
        a = 1;
    else
        a = 0;
    end
    
    mn = [mn a];
end

disp(' Binary information at Receiver :');
disp(mn);

% ========== Digital Signal Representation ==========
bit = [];
for n = 1:length(mn)
    if mn(n) == 1
        se = 5*ones(1,100);
    else
        se = zeros(1,100);
    end
    bit = [bit se];
end

t5 = bp/100:bp/100:100*length(mn)*(bp/100);
subplot(4,1,4);
plot(t5,bit,'LineWidth',2.5); grid on;
axis([0 bp*length(mn) -0.5 6]);
ylabel('amplitude(volt)');
xlabel('time(sec)');
title('Demodulated signal at receiver');
