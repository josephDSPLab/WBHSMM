function onset_x = OnsetsDetection(x,p0,hop,frame)
    onset_x = zeros(size(x));
    for pitch_i = 1:length(p0)
        if p0(pitch_i) == 0
            continue;
        end
        x_now = x((pitch_i-1)*hop+1:(pitch_i-1)*hop+frame);
        win_size = ceil(p0(pitch_i)/10);
        x_p = padarray(x_now,ceil(win_size/2),mean(x_now));
        onset=zeros(size(x_now));
        for i = 1:length(x_now)
            [m1,st] = min(x_p(i:i+win_size));
            m2 = max(x_p(i+st-1:i+win_size));
            onset(i) = m2-m1;
        end
        onset_candidate = diff(onset);
        [~,p] = sort(onset_candidate);
        peak = p(end-10:end);
        [~,onset_p] = max(abs(peak - length(x_now)/2));
        onset_p = peak(onset_p);
        onset(:) = 0;
        onset(onset_p) = 1;
        search_c = onset_p;
        if isempty(search_c)
            continue
        end
        while (search_c - p0(pitch_i) - ceil(win_size/2)) > 1
            search_c = search_c - p0(pitch_i);
            [~,p] = min(onset_candidate(search_c-ceil(win_size/2):search_c+ceil(win_size/2)));
            if isempty(p)
                continue
            end
            onset(search_c-ceil(win_size/2)+p+1) = 1;
        end
        search_c = onset_p;
        while (search_c + p0(pitch_i) + ceil(win_size/2)) < length(onset)
            search_c = search_c + p0(pitch_i);
            [~,p] = min(onset_candidate(search_c-ceil(win_size/2):search_c+ceil(win_size/2)));
            onset(search_c-ceil(win_size/2)+p+1) = 1;
        end
        onset_x((pitch_i-1)*hop+1:(pitch_i-1)*hop+frame) = onset;
%         plot(onset)
%         hold on 
%         plot(x_now)
%         hold off
%         k = waitforbuttonpress;
    end
    peak_candid_all = find(onset_x>0);
    for ii = 2:sum(onset_x>0)
        t_c = round(peak_candid_all(ii)/hop)+1;
        if peak_candid_all(ii) - peak_candid_all(ii-1) <median( p0(max(t_c-5,1):min(t_c+5,end)))/2
            onset_x(peak_candid_all(ii)) = 0;
        elseif peak_candid_all(ii) - peak_candid_all(ii-1) > median(p0(max(t_c-5,1):min(t_c+5,end)))*2
            onset_x(round((peak_candid_all(ii)+peak_candid_all(ii-1))/2)) = 1;
        end
    end

end