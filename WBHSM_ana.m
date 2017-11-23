function [para,onset] = WBHSM_ana(x,fs,segtime,hop,f0)
    %block-wise process
    fbins = 8192;
    framesize = round(segtime*fs);
    win = hamming(framesize);

    xblock = zeros([framesize,floor((length(x)-framesize)/hop)+2]);
    p0 = round(fs./f0);
    p0(p0 == inf) = 0;
    onset = OnsetsDetection(x,p0,hop,framesize);
    mode = 0;
    if mode ==0  % Not perform well
        [S_n,p_info] = perdiz(x,p0,onset,hop);
        para.S_n = S_n;
        para.p_info = p_info;
        para.mode = 'periodize';
        xblock(:,end) = [x((size(xblock,2)-1)*hop+1:end);zeros([framesize - length(x((size(xblock,2)-1)*hop+1:end)),1])].*win;
    else  % upsampling : failed currently
        [S_n,p_info] = ups(x,p0,onset,hop); 
        para.S_n = S_n;
        para.p_info = p_info;
        para.mode = 'ups';
        xblock(:,end) = [x((size(xblock,2)-1)*hop+1:end);zeros([framesize - length(x((size(xblock,2)-1)*hop+1:end)),1])].*win;
    end

end
function [S_n,p_info] = ups(x,p0_t,onset,hop)
    peak_candid_all = find(onset>0);
    p0_sup = find(p0_t > 0);
    p_info = zeros([2,sum(onset > 0)]);
    for ii = 1:sum(onset>0)
        onset_p = peak_candid_all(ii);
        t_c = round(peak_candid_all(ii)/hop)+1;
        [~,t_pos] = min(abs(p0_sup - t_c));
        p0 = p0_t(p0_sup(t_pos));
        p_concat = round(p0*1);
        p_info(1,ii) = p0;
        oneperiod = x(onset_p-round(p_concat/2):onset_p-round(p_concat/2)+p_concat-1);

        %%%
        fftbins = 8192;
        fmax = 6000;
        faxes = linspace(0,44100*fftbins/length(oneperiod),fftbins);
        maxbinspos = sum(faxes<fmax);
        fmax = faxes(maxbinspos);
        %%%
        win_ana = hamming(fftbins);
        x_r_win = resample(oneperiod,fftbins,length(oneperiod)).*win_ana;
        plot(x_r_win)
        dft_xr = fft(x_r_win,fftbins);
        p_info(2,ii) = fftbins;
        %dft_xr = dft_xr(1:fftbins/2+1);
        Xm = abs(dft_xr(1:maxbinspos));
        f0 = 44100/p0;
        ParaNum = floor((fmax-f0)/f0);  
        [a_k,f_k] = HarmonPeakFind(Xm,f0,22050*fftbins/length(oneperiod)/(fftbins/2),faxes,ParaNum);
        Xtheta = angle(dft_xr(1:maxbinspos));
        theta_k = Xtheta(f_k);
        S_n{ii} = {a_k,f_k,theta_k};
    end
end
function [S_n,p_info] = perdiz(x,p0_t,onset,hop)

    peak_candid_all = find(onset>0);
    p0_sup = find(p0_t > 0);
    p_info = zeros([2,sum(onset > 0)]);
    for ii = 1:sum(onset>0)

        onset_p = peak_candid_all(ii);
        t_c = round(peak_candid_all(ii)/hop)+1;
        [~,t_pos] = min(abs(p0_sup - t_c));
        p0 = p0_t(p0_sup(t_pos));
        num_period = 7;%ceil(8192/p0*5);
        num_period = num_period + (1-mod(num_period,2));
        p_concat = round(p0*1.5);
        p_info(1,ii) = p_concat;
        win = tukeywin(p_concat,0.25);
        oneperiod = x(onset_p-round(p_concat/2):onset_p-round(p_concat/2)+p_concat-1).*win;
        x_periodized = zeros([round(p0*(num_period+0.5)),1]);
        %calculate anchor point
        anc = round([1:p0:p0*(num_period-1)+1]);
        for pulse_pos = anc
            x_periodized(pulse_pos:pulse_pos+p_concat-1) = x_periodized(pulse_pos:pulse_pos+p_concat-1) + oneperiod;
        end

        win_ana = hamming(anc(end)-anc(1)+1);
        p_info(2,ii) = anc(end)-anc(1)+1;
        x_r_win = x_periodized(anc(1)+round(p_concat/2):anc(end)+round(p_concat/2)).*win_ana;
        %%%
        fftbins = length(x_r_win);
        fmax = 6000;
        faxes = linspace(0,44100,fftbins);
        maxbinspos = sum(faxes<fmax);
        fmax = faxes(maxbinspos);
        %%%
        dft_xr = fft(x_r_win,fftbins);
        %dft_xr = dft_xr(1:fftbins/2+1);
        Xm = abs(dft_xr(1:maxbinspos));
        f0 = 44100/p0;
        ParaNum = floor((fmax-f0)/f0);  
        [a_k,f_k] = HarmonPeakFind(Xm,f0,22050/(fftbins/2),faxes,ParaNum);
        Xtheta = angle(dft_xr(1:maxbinspos));
        theta_k = Xtheta(f_k);
        S_n{ii} = {a_k,f_k,theta_k};
    end
end
function [a_k,f_k] = HarmonPeakFind(Xm,f0,fqres,faxes,ParaNum)
    winsize = round((f0*0.1)/fqres);
    a_k = abs(Xm(1));
    f_k = 1;
    for i = 1:ParaNum
        fi_pos = sum(f0*i > faxes);
        search_fbins = Xm(fi_pos-winsize:fi_pos+winsize);
        [apeak,fpeak] = max(search_fbins);
        a_k = [a_k;search_fbins(search_fbins>apeak*0.9)];
        f_k = [f_k;find(search_fbins>apeak*0.9)+fi_pos-winsize-1];
    end
end