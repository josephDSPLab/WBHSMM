function x_syn = WBHSM_syn(x,fs,para,onset)
    S_n = para.S_n;
    p_info = para.p_info;
    g_idx = [-10:10]';
    Ga = exp(-g_idx.^2);%fftshift(abs(fft(hamming(8192))));
    anglewin = ones(size(Ga));
    x_syn = zeros(size(onset));
    onset_p = find(onset > 0);
    x_syn(1:onset_p(1)) = x(1:onset_p(1));
    x_syn(onset_p(end):end) = x(onset_p(end):end);
    for ii = 1:length(S_n)
        p0 = p_info(1,ii);
        a_k = S_n{ii}{1};
        f_k = S_n{ii}{2};
        theta_k = S_n{ii}{3};
        k = 0.25;
        k2 = 1+k;
        
        fftbuffer = zeros([p_info(2,ii),1]);
        fftbuffer(f_k) = a_k;
        fftbuffer(end-f_k(2:end)+2) = a_k(2:end);
        fftbuffer = conv(fftshift(fftbuffer),Ga,'same');
        anglebuffer = zeros(size(fftbuffer));
        anglebuffer(f_k) = theta_k;
        anglebuffer(end-f_k(2:end)+2) = -theta_k(2:end);
        anglebuffer = conv(fftshift(anglebuffer),anglewin,'same');
        fftbuffer = fftbuffer.*exp(anglebuffer.*1j);
        x_periodized = real(ifft(fftshift(fftbuffer)));
        x_periodized = fftshift(x_periodized);
        if strcmp(para.mode,'periodize')
            win_syn = tukeywin(round(p0*k2),k/2);
            st = round(p_info(2,ii)/2)-round(p0*k2/2);
            oenperiod = x_periodized(st:st + round(p0*k2 ) - 1).*win_syn;
            x_syn(onset_p(ii):onset_p(ii)+round(p0*k2 )-1) =  x_syn(onset_p(ii):onset_p(ii)+round(p0*k2 )-1)+oenperiod;
        else
            win_syn = tukeywin(p0,k/2);
            oenperiod = resample(x_periodized,p0,p_info(2,ii)).*win_syn;
            x_syn(onset_p(ii):onset_p(ii)+p0-1) =  x_syn(onset_p(ii):onset_p(ii)+p0-1)+oenperiod;
        end
        
    end
end