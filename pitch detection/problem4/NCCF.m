function [F0,m] = NCCF(x,fs,segtime,hop,fmin,fmax)
%Homework 2
%   Joseph Liao
% NCCF F0 contour
lowestp = round(fs/fmin);
pmax = round(fs/fmax);
% Extract frequency information
framesize = round(segtime*fs);

xblock = zeros([framesize,floor((length(x)-framesize)/hop)+2]);
for i = 1:size(xblock,2)-1
    xblock(:,i) = x((i-1)*hop+1:(i-1)*hop+framesize);
end
xblock(:,end) = [x((size(xblock,2)-1)*hop+1:end);zeros([framesize - length(x((size(xblock,2)-1)*hop+1:end)),1])];

nccf = zeros([framesize,floor((length(x)-framesize)/hop)+2]);
for j = pmax:lowestp
    nccf(j,:) = sum(xblock(1:lowestp,:).*xblock(j:j+lowestp-1,:),1);
    nccf(j,:) = nccf(j,:)./( sqrt(sum(xblock(1:lowestp,:).*xblock(1:lowestp,:))).*sqrt(sum(xblock(j:j+lowestp-1,:).*xblock(j:j+lowestp-1,:))) );
end
nccf(1:pmax,:) = 0;
%nccf = nccf./repmat(max(nccf),[size(nccf,1),1]);
[merits,p]= sort(abs(nccf),'descend');
p = p(1:3,:);
F0 = fs./p;
m = merits(1:3,:);
end