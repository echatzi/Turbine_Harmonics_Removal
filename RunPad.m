
clear; close all; clc

%% Load data
load('C:\Users\chatz\Downloads\6200135_2024-02-07.mat')
TD_resampled = addMBCChannelsToTD(TD_resampled);

idx = 4e5:6e5;

fs = 10;  % Sampling frequency [Hz]

x = TD_resampled.BLSSMV_BldC_SS.value(idx,end);
x = x(:);
x = x - mean(x);

% Rotor 1P frequency [Hz]
% Assumes generator speed is in rpm
gearRatio = 126.793;
f1P = TD_resampled.GenSdSpd.value(idx,end) ...
      / (gearRatio*60);
f1P = f1P(:);

f1Pmean = mean(f1P);

fprintf('Mean 1P frequency: %.4f Hz\n',f1Pmean);
fprintf('1P standard deviation: %.4f Hz\n',std(f1P));

%% Optional 1P notch filtering
remove1P = false;  % Change to true to remove 1P

xFiltered = x;

if remove1P
    % Fixed notch, appropriate for nearly constant speed
    notchBW = max(0.03, 0.20*f1Pmean); % Hz

    w0 = f1Pmean/(fs/2);
    bw = notchBW/(fs/2);

    [b,a] = iirnotch(w0,bw);

    xFiltered = filtfilt(b,a,x);
end

%% STFT parameters
windowDuration = 60; % seconds
windowLength = round(windowDuration*fs);
overlapLength = round(0.90*windowLength);

fftLength = 2048;
nColors = 128;

%% Compute STFT
[S,F,T] = stft(xFiltered,fs, ...
    Window=hann(windowLength,"periodic"), ...
    OverlapLength=overlapLength, ...
    FFTLength=fftLength, ...
    FrequencyRange="onesided");

%% Relative magnitude in dB
Sdb = 20*log10(abs(S)+eps);
Sdb = Sdb-max(Sdb(:));

dynamicRange = [-80 0];

%% 2D spectrogram
figure

imagesc(T,F,Sdb)
axis xy

colormap(turbo(nColors))
clim(dynamicRange)
ylim([0 5])

colorbar
xlabel('Time [s]')
ylabel('Frequency [Hz]')
title('STFT - 2D')

% Overlay the instantaneous 1P frequency
tRaw = (0:length(x)-1)'/fs;
hold on
k = 1:100:length(x);
hold on

orders = [1 2 3 6 9];
colors = lines(length(orders));

for i = 1:length(orders)
    plot(tRaw(k), orders(i)*f1P(k), ...
        'Color', colors(i,:), ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('%dP',orders(i)));
end

legend('show')
hold off

%% 3D spectrogram
% Reduce plotting density to improve graphics performance
stepT = max(1,ceil(length(T)/400));
idxT = 1:stepT:length(T);

idxF = find(F<=5);
idxF = idxF(1:3:end);

figure

surf(T(idxT),F(idxF), ...
     Sdb(idxF,idxT), ...
     'EdgeColor','none');

colormap(turbo(nColors))
clim(dynamicRange)

xlabel('Time [s]')
ylabel('Frequency [Hz]')
zlabel('Relative magnitude [dB]')
title('STFT - 3D')

ylim([0 5])
zlim(dynamicRange)

view(40,45)
grid on
colorbar
