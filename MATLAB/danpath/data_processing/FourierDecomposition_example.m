%FourierDecomposition_example.m

clear;close('all');clc

%% Generate signal with two frequencies (50, 120 Hz) + noise, sampled at 1000 Hz
Fs = 1000;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = 1000;                     % Length of signal
t = (0:L-1)*T;                % Time vector
% Sum of a 50 Hz sinusoid and a 120 Hz sinusoid
freq1 = 50; %[hZ]
amp1 = .7;  %[-]
freq2 = 120; %[hZ]
amp2 = 1;  %[-]
x = amp1*sin(2*pi*freq1*t) + amp2*sin(2*pi*freq2*t); 
y = x + 2*randn(size(t));     % Sinusoids plus noise

figure(1)
plot(Fs*t(1:50),y(1:50))
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('time (milliseconds)')

%% Calculate power spectrum using fft
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
figure(2)
plot(f,2*abs(Y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')