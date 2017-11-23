function [f0_detection,C,back_trace] = DP_F0track(F0,F0_shc,m,m_shc,nlfer)
%Homework 2
%   Joseph Liao
% Dynamic programming find best F0 track
C = 1-m;   %initialize cost value with each time stamp
meanf0_1 = sum(F0(:,1))/sum(F0(:,1)>0);
meanf0_2 = sum(F0_shc(:,1))/sum(F0_shc(:,1)>0);
back_trace = inf(size(m));
for d = 2:size(F0,1)
        f1 = repmat(F0(d-1,:)',[1 size(F0,2)]);
        f2 = repmat(F0(d,:),[size(F0,2) 1]);
        tr_cost = transcost(f1,f2,meanf0_1,abs(nlfer(d)-nlfer(d-1)));
        [min_c,back_trace(d,:)] = min(0.15*tr_cost.*repmat(C(d-1,:)',[1 size(F0,2)]));
        f1 = repmat(F0(d,:)',[1 size(F0,2)]);
        f2 = repmat(F0_shc(d,:),[size(F0,2) 1]);
        tr_cost = transcost(f1,f2,meanf0_2,abs(nlfer(d)-nlfer(d-1)));
        min_c2= min(tr_cost);
        C(d,:) = C(d,:)+min_c2 + min_c;
end
% Backtrace routine
[back_cost, p] = min(C(end,:));
f0_detection = F0(end,p);
for b_index = size(C,1):-1:2
    p = back_trace(b_index,p);
    f0_detection = [F0(b_index-1,p),f0_detection];
end
end
function cost = transcost(f1,f2,stdf0,diff_nlfer)
    cost = 0.9*abs(f1-f2)./stdf0;
    cost(f1==0 | f2==0) = 0.5*min(1,diff_nlfer);
    cost(f1==0 & f2 ==0) = 0.1;
end