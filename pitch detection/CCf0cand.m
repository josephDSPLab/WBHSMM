function [F0,t,m] = CCf0cand(x,fs,segtime,hop,fmin,fmax)
%Homework 2
%   Joseph Liao
% Cesptral F0 contour

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
Sm = abs(S);     % preserve the magnitude value
% show the interesting region

CC = ifft(log10(Sm));
CC = CC(1:floor(fbins/2),:);
%liftering cut-off qrefrency
qc = 30;
pmax = round(fs/fmax);
pmin = round(fs/fmin);
CCt = CC;
CCt(qc:end,:) = 0;
CCp = CC;
CCp(1:pmax,:) = 0;
CCp(pmin:end,:) = 0;
CCpm = abs(CCp);
%refinement
CCpm = medfilt2(CCpm,[3 1]);
%CCpm = CCpm./repmat(max(CCpm),[size(CCpm,1),1]);
[merits,p] = sort(CCpm,'descend');
p = p(1:3,:);
%CCpm(CCpm<repmat(maxv*0.25,[size(CCpm,1),1])) = 0;
%peak_num = sum(CCpm>0,1);
F0 = fs./p;

m = merits(1:3,:);
end
