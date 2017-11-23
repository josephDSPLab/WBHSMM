function [F0,m,t,uv] = SHCf0cand(x,fs,segtime,hop,fmin,fmax)
%Homework 2
%   Joseph Liao
% Spectral Harmonic Correlation F0 contour

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
NH = 3;
f_axes = linspace(0,fs/2,fbins/2+1);
pmin = sum(f_axes < fmin)+1;
pmax = sum(f_axes < fmax)+1;
SHC = zeros(size(Sm));
for jj = pmin:pmax
    %calculate window range
    f_delta = 40;
    for har_c = 1:NH
        f_selected = har_c * f_axes(jj);
        f_win = abs(f_axes-f_selected) < f_delta/2;
        if isempty(f_win)
            continue
        end
        SHC(jj,:) = SHC(jj,:) + sum(Sm(f_win,:));
    end
end
norm_c = max(SHC);
SHC = SHC./repmat(norm_c,[size(SHC,1) 1]);
[merits,p] = sort(SHC,'descend');
uv = merits(1,:) < 0.2;

p = p(1:3,:);
F0 = f_axes(p);
m = merits(1:3,:);
end
