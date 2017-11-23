function nlfer = NLFER(x,fs,segtime,hop,fmin,fmax)
%Homework 2
%   Joseph Liao
% NLFER

% Extract frequency information
fbins = 8192;
framesize = round(segtime*fs);
win = hamming(framesize);
xblock = zeros([framesize,floor((length(x)-framesize)/hop)+2]);
for i = 1:size(xblock,2)-1
    xblock(:,i) = x((i-1)*hop+1:(i-1)*hop+framesize).*win;
end
xblock(:,end) = [x((size(xblock,2)-1)*hop+1:end);zeros([framesize - length(x((size(xblock,2)-1)*hop+1:end)),1])].*win;
t = [0:size(xblock,2)-1]*hop;
S = fft(xblock,fbins);
Sm = abs(S(1:fbins/2+1,:));     % preserve the magnitude value
f_axes = linspace(0,fs/2,fbins/2+1);
pmin = sum(f_axes < fmin*2)+1;
pmax = sum(f_axes < fmax)+1;
nlfer = sum(Sm(pmin:pmax,:),1);
nlfer = nlfer / mean(nlfer);
nlfer(nlfer < 0.75) = 0;
end