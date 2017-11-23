close all;
[x, fs] = audioread('mic_F01_sa1.wav');
Isplotcand = 1;
gt = textread('ref_F01_sa1.f0');
gt = gt(:,1);
f0_detection = SpecTempF0Track(x,fs,Isplotcand);
x1 = linspace(0,length(x)/fs,length(x));
x2 = linspace(0,length(x)/fs,length(f0_detection));
x3 = linspace(0,length(x)/fs,length(gt));

figure
plotyy(x1,x,[x2',x2'],[f0_detection',interp1(x3,gt,x2)'])
legend('Signal','My result','Ground truth')
set(gca,'FontSize',16)
xlabel('Time (s)','fontsize',20)