%Homework 2 MFCCs
%   Joseph Liao
% Using triangular filter bank to calculate MFCC

% load data
[x, fs] = audioread('mic_F01_sa1.wav');

% Pre-emphasis
alpha = 0.95;
a = 1;
b = [1 -alpha];
s = x;%s = filter(b,a,x);

% Extract frequency information
fbins = 1024;
win = hamming(0.05*fs);
overlap = round(0.03*fs);
[S,f_axis,t] = spectrogram(s,win,overlap,fbins,fs);
Sm = abs(S);     % preserve the magnitude value
f_axis = linspace(0,fs/2,fbins/2+1);

% Set Mel frequency parameter
Melbins = 40;
fmin = 100;
fmax = 4000;

%calculate center point of each triangle filter bank
h2mel = @(f) 1000/log10(2).*log10(1+f./1000);
mel2f = @(f) 1000.*(10.^(f.*log10(2)/1000) - 1);
        % %simple test
        % melmin = h2mel(fmin);
        % fmin2 = mel2f(melmin);
tribank = linspace(h2mel(fmin),h2mel(fmax),Melbins+2);
fc = mel2f(tribank);     %convert center point back into Hertz scale

%initialize filterbank matrix (size : Melbins x fbins)
Htri = zeros([Melbins,length(f_axis)]);
for i = 1:Melbins
    Htri(i,:) = Htri(i,:) + (f_axis>=fc(i)&f_axis<=fc(i+1)).*(f_axis-fc(i))/(fc(i+1)-fc(i));
    Htri(i,:) = Htri(i,:) + (f_axis>=fc(i+1)&f_axis<=fc(i+2)).*(fc(i+2)-f_axis)/(fc(i+2)-fc(i+1));
end
%plot tri- filter bank
        % figure
        % hold on
        % for i = 1:Melbins
        %     plot(Htri(i,:),'b')
        % end
        % hold off
%
Smel = Htri* Sm;
dctm = @( L, Melbins )(cos( repmat([0:L-1].',[1,Melbins]).* repmat(pi*([1:Melbins]-0.5)/Melbins,[L,1]) ) );

DCT = dctm( 20, Melbins );
MFCCs = DCT * log10(Smel);
imagesc( [1:size(MFCCs,2)], [0:size(MFCCs,1)-1], abs(MFCCs) );
axis( 'xy' );