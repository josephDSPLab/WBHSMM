function f0_detection = SpecTempF0Track(x,fs,Isplotcand)
    F = [50 100 1450 1550];  % band limits
    A = [0 1 0];                % band type: 0='stop', 1='pass'
    dev = [0.0001 10^(0.1/20)-1 0.0001]; % ripple/attenuation spec
    [M,Wn,beta,typ] = kaiserord(F,A,dev,fs);  % window parameters
    b = fir1(M,Wn,typ,kaiser(M+1,beta),'noscale'); % filter design
    x = filter(b,1,x);

    n_candidate = 3;
    segtime = 0.032;
    fmin = 100;
    fmax = 1000;
    hop = round(0.01*fs);
    x_nonlinear = abs(x);
    nlfer = NLFER(x_nonlinear,fs,segtime,hop,fmin,fmax);
    [F0_nccf,m_nccf] = NCCF(x_nonlinear,fs,segtime,hop,fmin,fmax);
    [F0_shc,m_shc,t,uv] = SHCf0cand(x_nonlinear,fs,segtime,hop,fmin,fmax);
    voice_seg = nlfer>0 & uv == 0;
    F0_c5 = F0_nccf.*repmat(voice_seg,[n_candidate,1]);
    F0_candidate = [F0_c5'];
    F0_candidate = medfilt1(F0_candidate,5,1);
    F0_shc = medfilt1(F0_shc,7,[],2);
    m_candidate = [m_nccf'];
    x_ax = repmat(linspace(0,length(x)/fs,size(F0_candidate,1)),[n_candidate,1]);
    if Isplotcand
        plot(x_ax',F0_candidate,'o')
        legend('Candidate 1','Candidate 2','Candidate 3')
        set(gca,'FontSize',16)
        xlabel('Time (s)','fontsize',20)
        ylabel('Frequency (Hz)','fontsize',20)
    end
    [f0_detection,C,back] = DP_F0track(F0_candidate,F0_shc',m_candidate,m_shc',nlfer);
end