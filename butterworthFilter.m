%{
    Author: Sarker Nadir Afridi Azmi
%}

% |Ha(jOmega)|^2 = 1 / (1 + (jOmega/jOmegac)^2n

% Clear the workspace
clear;

% 1. Audio Frequency Analysis
% ---------------------------
[sampledData, Fs] = audioread("noisyaudio.wav");
dft = fft(sampledData);
% Find the correct dimensions
dft = dft';
% (upper limit - lower limit) / total sample size
f = (-Fs/2):(Fs)/(length(sampledData) - 1):(Fs/2);

magDftw0 = fftshift(abs(dft));
subplot(2, 2, 1);
plot(f,magDftw0);
title("DFT Magnitude vs Frequency");
xlabel("Frequency [Hz]"), ylabel("Magnitude");

% Normalized dft - kinda confused right now.
normalizedLogPlotDft = 20.*log10(abs(dft)/max(abs(dft)));
subplot(2, 2, 2);
plot(f, fftshift(normalizedLogPlotDft));
title("Normalized log plot of Magnitude vs Frequency");
xlabel("Frequency [Hz]"), ylabel("Magnitude [dB]");

% 2. Filter Design
% ----------------
wp = mapTpW(1073, Fs/2);            % pass band
ws = mapTpW(2267, Fs/2);            % stop band
delp = -1;                          % 1 - delp (pass band)
dels = -80;                         % stop band

% Equations
kp = (1/(convertFromDB(delp))^2) - 1;
ks = (1/convertFromDB(dels)^2) - 1;

% Find the order of the butterworth filter
N = (0.5*log10(ks/kp))/log10(transformBilinearly(ws)/transformBilinearly(wp));
N = ceil(N);

% There will be no aliasing
% So, we can meet the stop band requirement
% Calculate corrected cut off frequency for the analog butterworth filter
BigOmegaC = transformBilinearly(ws)/(ks^(1/(2*N)));

% Equivalent cut-off frequency
Fc = atan(BigOmegaC/2);

% Logarithmic gain of frequency response of analog filter
HaS = 20*log10(sqrt(1./(1 + (f/Fc).^(2*N))));
subplot(2, 2, 3);
plot(f,fftshift(abs(HaS)));
title("Logarithmic gain of frequency response");
xlabel("Frequency [Hz]"), ylabel("||Ha(s)|");

% 3. Filter implementation
% ------------------------
%{
    This might seem weird, but it pertains to what is said in 3a
    Fc is in radians/sec and so, the span of the digital filter is [-pi,
    pi] instead of [-Fs/2, Fs/2]
%}
Wn = Fc/(pi);

[b,a] = butter(N, Wn);

filteredAudio = filter(b,a,sampledData);
filteredAudioDft = fft(filteredAudio);
subplot(2, 2, 4);
plot(f,fftshift(abs(filteredAudioDft)));
title("Filtered audio DFT Magnitude vs Frequency");
xlabel("Frequency [Hz]"), ylabel("Magnitude");

sound(filteredAudio, Fs);
audiowrite("filteredAudio.wav", filteredAudio, Fs);

% Maps the frequency to angular frequency
function res = mapTpW(f, Fs)
    res = (f/Fs)*pi;
end

% Finds the bilinear transformation
function bOmega = transformBilinearly(sOmega)
    bOmega = 2*tan(sOmega);
end

function x = convertFromDB(xdB)
    x = 10^(xdB/20);
end