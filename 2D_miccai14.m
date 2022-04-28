%% ---------- Mohammad Arafat Hussain ---------------------%%
% comments and concerns to m.arafathussain@gmail.com
%
% ------------------------------------------------------------
clc;
clear all;
close all; 

getd('C:\Users\arafat\Dropbox\AM2D'); % change directory accordingly

volume1 = loadvol('W:\data\US\Elastography Data\feb 11 2015\05sin.vol',33);
volume2 = loadvol('W:\data\US\Elastography Data\feb 11 2015\05sin.vol',35);

rf1 = volume1; clear volume1;
rf2 = volume2; clear volume2;

%%
% B-mode
% env = envelope(rf1);
% option = 'mu';
% fac = 255;
% B = LogComp(env,option,fac);

% imagesc(B);colormap gray
%%
rf1 = abs(hilbert(rf1));
rf2 = abs(hilbert(rf2));

[x,y] = size(rf1);
rf = rf1.^2;
ttb = zeros(x,y);
ttb(1,:) = rf(1,:);

for i = 2:x
    ttb(i,:) = rf(i,:) + ttb(i-1,:);
    if i == x
        for j = 1:y
         ttb(:,j) = ttb(:,j)./max(ttb(:,j));
        end
    end
end
rf1 = rf.^ttb;

rf = rf2.^2;
ttb = zeros(x,y);
ttb(1,:) = rf(1,:);

for i = 2:x
    ttb(i,:) = rf(i,:) + ttb(i-1,:);
    if i == x
        for j = 1:y
         ttb(:,j) = ttb(:,j)./max(ttb(:,j));
        end
    end
end
rf2 = rf.^ttb;

%% strain estimation
IRF = [-20 0];
alfa_DP = 0.15;
midA = 64;
alfa = 1;  
beta = 1;
gamma = 0.005;
T = 0.4;
f_0 = 5;
f_s = 40;
wDIff = 5;

str = am2d_expph(rf1,rf2,IRF,alfa_DP,midA,alfa,beta,gamma,T,f_0,f_s,wDIff);
% imagesc(str);colormap gray

%% strain processing
str = (-1).*str;
minI = median(str(:));
str = (abs(str + abs(minI))).^1;
maxIm = max(str(:));
str = str/maxIm;
str = imadjust(str);
sn = str;

%% RF processing
bm = abs(hilbert(rf1));
rf = bm.^2;
rf = medfilt2(rf);
maxIm = max(rf(:));
rf = rf/maxIm;
rf = imadjust(rf);
[aa,bb] = size(rf);
rn = rf(43:size(rf,1)-43,11:size(rf,2)-10);
% figure;imagesc(B(43:size(rf,1)-43,11:size(rf,2)-10));colormap gray;axis square;hold on
% clear rf;

% Fused Map
final = 0.5*rn + 0.5*sn;

t = zeros(1,size(final,2));
tmp = max(final);

for i = 1:size(final,2)
     t(i) = max(find(tmp(i) == final(:,i)));
end

a = 20;

for l = 1:size(final,2)-a+1
    win_to_fit = t(l:l+a-1);
    len = length(win_to_fit);
    linearCoef = polyfit((1:len)',win_to_fit',1);
    t(l:l+a-1) = polyval(linearCoef,(1:len)');           
end
t = round(t.*(364/2048));   % Final 2D bone
% plot(t)
t = fliplr(t);