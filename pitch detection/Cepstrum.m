%Homework 2 MFCCs
%   Joseph Liao
% Using triangular filter bank to calculate MFCC
close all
% load data
[x, fs] = audioread('sample.wav');


% Extract frequency information
fbins = 1024;
framesize = round(0.02*fs);
win = hamming(framesize);
over_r = 0.75;
hop = round((1-over_r)*framesize);
xblock = zeros([framesize,ceil((length(s)-framesize)/overlap)]);
for i = 1:size(xblock,2)-1
    xblock(:,i) = x((i-1)*hop+1:(i-1)*hop+framesize).*win;
end
xblock(:,end) = [x((i-1)*framesize+1:end);zeros([framesize - length(x((i-1)*framesize+1:end)),1])].*win;
S = fft(xblock,fbins);
Sm = abs(S);     % preserve the magnitude value
% show the interesting region

CC = ifft(log10(abs(S)));
CC = CC(1:floor(fbins/2)+(1-mod(fbins,2)),:);
%liftering cut-off qrefrency
qc = 30;
CCt = CC;
CCt(qc:end,:) = 0;
CCp = CC;
CCp(1:qc,:) = 0;
plot(abs(S(1:floor(fbins/2)+(1-mod(fbins,2)),size(S,2)/2)))
figure
plot(abs(CCp(:,size(S,2)/2)))
[~,p] = max(abs(CCp));