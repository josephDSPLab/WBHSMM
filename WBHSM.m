function [para,re_x] = WBHSM_ana(x,fs,segtime,f0)
    %STFT
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
    %Periodization
end
function x_per = perdiz(x,f0)
    
end