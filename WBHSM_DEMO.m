 %WBHSM
close all;
addpath('pitch detection');
[x, fs] = audioread('femalesingingjazz.wav');
Isplotcand = 0;
gt = textread('ref_F01_sa1.f0');
gt = gt(:,1);
f0_detection = SpecTempF0Track(x,fs,Isplotcand);
f0_detection(f0_detection == 0) = round(mean(f0_detection(f0_detection>0)));
f0_detection([1 end]) = 0;

[para,onset] = WBHSM_ana(x,fs,0.032,0.01*fs,f0_detection);
x_syn = WBHSM_syn(x,fs,para,onset);