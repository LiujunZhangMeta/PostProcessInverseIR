clc;
clear all;
close all;
% Default regularization settings
config.f1 = 200;      % Low cut-off of the headphones reproduction range
config.f2 = 16000;    % High cut-off of the headphones reproduction range
config.srate = 48000; % Sampling rate
config.taps = 2048;   % Number of frequency taps of the inverse filter
config.order1 = 5;    % Order of the high-pass Butterworth filter (def = 5)
config.order2 = 5;    % Order of the low-pass Butterworth filter (def = 5)
config.option = 3;    % 0: do not preserve low-frequencies
config.max_amp = 20;  % Maximum amplification allowed by the regularization
config.k = 0.5;       % Smoothing window length for notch regularization
config.phi = 0;       % Factor used to alter the phase of the
config.targetIR = 1;  % Target EQ curve (default = 1, which results in flat EQ)
[IR,fs] = audioread('headphoneIR4144-5.wav');
invIR = autoreg(IR,config);
%% plot the inverse IR and the original input IR
figure(2)
subplot(2,1,1)
plot(IR(:,1),LineWidth=1);hold on;plot(IR(:,2)+0.8,LineWidth=1);
xlim([0,2048]);%ylim([-40,15]);
grid on; xlabel('Sample');
set(gcf,'color','w');ylabel('Amplitude');hold on
title('Impulse response');
subplot(2,1,2)
plot(invIR(:,1),LineWidth=1);hold on;plot(invIR(:,2)+0.8,LineWidth=1);
xlim([0,2048]);%ylim([-40,15]);
grid on; xlabel('Sample');
set(gcf,'color','w');ylabel('Amplitude');hold on
title('Inverse Impulse Repsonse');




