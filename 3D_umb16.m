%% ---------- Mohammad Arafat Hussain ---------------------%%
% comments and concerns to m.arafathussain@gmail.com
%
% ------------------------------------------------------------
clc;
clear all;
close all; 

getd('C:\Users\arafat\Dropbox\ZPclustering'); % change directory accordingly
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

% 3D RF loading
volumeRF = loadvol('W:\data\US\Elastography Data\feb 11 2015\04mul.vol',2);
RF = double(volumeRF(43:aa-43,11:bb-10,38:64));

% Surface Growing
volume = loadvol('W:\data\US\Elastography Data\feb 11 2015\05mulB.vol',1);
USenv = double(volume(:,78:3:406,38:64));
% save_untouch_nii(make_ana(USenv,[.3 .3 .3], [0 0 0]),'Invivo2-Bmode');

[aa,bb,cc] = size(USenv);
rf = USenv.^2;
ttb = zeros(aa,bb,cc);
ttb(1,:,:) = rf(1,:,:);

for i = 2:aa
    ttb(i,:,:) = rf(i,:,:) + ttb(i-1,:,:);
    if i == aa
        for j = 1:bb
            for k = 1:cc
                ttb(:,j,k) = ttb(:,j,k)./max(ttb(:,j,k));
            end
        end
    end
end
USenv = USenv.^ttb;


USenv = USenv(:,1:108,:);
USenv = USenv./max(USenv(:));

fr = USenv.^1;%clear USenv;
[x,y,z] = size(fr);
idx = round(z/2);
ci = zeros(bb,cc);

span = 20;
Ecost = ones(2*span+1,y,z);
indx = zeros(2,y,z);
final_bone = zeros(y,z);
bs = zeros(z,1);
alfa = 2;
gamma = 0.8;

nny = 1;
nnz = 5;
ch = 1;
cv = 1;

for j = 1:y
    % clf
    xx = round(t(j));
    for k = idx:z
        if k == idx
            beta = 0;
        else
            beta = 0.5;
        end
        %c = ones(2*span+1,1);
        tmp1 = ones(2*span+1,1);
        tmp2 = ones(2*span+1,1);
        ind1 = xx - span; if ind1 < 1; ind1 = 1; end
        ind2 = xx + span; if ind2 > x; ind2 = x; end
        
        % intensity cost
        indx(1,j,k) = ind1;
        indx(2,j,k) = ind2;
        col = fr(:,j,k);
        st = mean(col(1:ind1));
        sb = mean(col(ind2:x));
        for i = ind1:ind2
            tmp1(i-ind1+1) = (1 - (fr(i,j,k) - st)^alfa - (fr(i,j,k) - sb)^alfa);
            tmp2(i-ind1+1) = (abs(xx - i)/span)^alfa;
            %c(i-ind1+1) = beta*(1 - (fr(i,j,k) - st)^alfa - (fr(i,j,k) - sb)^alfa) + (1 - beta)*0.5*(abs(xx - i)/span)^alfa;
        end
        tmp1 = tmp1./max(tmp1);
        [ci(j,k),~] = min(tmp1);
        tmp2 = tmp2./max(tmp2);
        c = tmp1 + beta*tmp2;clear tmp1;clear tmp2;
        [~,in] = min(c);
        xx = in + ind1 - 1;
        Ecost(:,j,k) = c;
        clear c;
    end
    
    xx = round(t(j));
    for k = idx-1:-1:1
        if k == idx-1
            beta = 0;
        else
            beta = 0.5;
        end
        %c = ones(2*span+1,1);
        tmp1 = ones(2*span+1,1);
        tmp2 = ones(2*span+1,1);
        ind1 = xx - span; if ind1 < 1; ind1 = 1;end
        ind2 = xx + span; if ind2 > x; ind2 = x;end
        
        % Intensity cost
        indx(1,j,k) = ind1;
        indx(2,j,k) = ind2;
        col = fr(:,j,k);
        st = mean(col(1:ind1));
        sb = mean(col(ind2:x));
        for i = ind1:ind2
            tmp1(i-ind1+1) = (1 - (fr(i,j,k) - st)^alfa - (fr(i,j,k) - sb)^alfa);
            tmp2(i-ind1+1) = (abs(xx - i)/span)^alfa;
            %c(i-ind1+1) = beta*(1 - (fr(i,j,k) - st)^alfa - (fr(i,j,k) - sb)^alfa) + (1 - beta)*0.5*(abs(xx - i)/span)^alfa;
        end
        tmp1 = tmp1./max(tmp1);
        [ci(j,k),~] = min(tmp1);
        tmp2 = tmp2./max(tmp2);
        c = tmp1 + beta*tmp2;clear tmp1;clear tmp2;
        [~,in] = min(c);
        xx = in + ind1 - 1;
        Ecost(:,j,k) = c;
        clear c;
    end  
end

for j = 1:y
    for k = 1:z
        
        j1 = j - nny; if j1 < 1; j1 = 1; end
        j2 = j + nny; if j2 > y; j2 = y; end
        k1 = k - nnz; if k1 < 1; k1 = 1; end
        k2 = k + nnz; if k2 > z; k2 = z; end
        
        w = zeros(j2-j1+1,k2-k1+1);
        ind_s = min(min(indx(1,j1:j2,k1:k2)));
        ind_e = max(max(indx(2,j1:j2,k1:k2)));
               
        tmp = ones(ind_e-ind_s+1,j2-j1+1,k2-k1+1);
         
        for  j11 = j1:j2
            for k11 = k1:k2
                if length(Ecost(:,j11,k11)) < 2*span+1
                    if indx(1,j11,k11) == 1
                        s = 1;
                        e = 2*span;
                    elseif indx(1,j11,k11) == x
                        s = (ind_e-ind_s+1) - 2*span + 1;
                        e = ind_e-ind_s+1;
                    end
                else
                    s = indx(1,j11,k11) - ind_s + 1;
                    e = s + 2*span;
                end
                tmp(s:e,j11-j1+1,k11-k1+1) = Ecost(:,j11,k11).*(exp(-abs(j11-j)*ch-abs(k11-k)*cv));
            end
        end
        
        cost = zeros(length(tmp(:,1,1)),1);
        
        for  j11 = j1:j2
            for k11 = k1:k2
                cost = cost + tmp(:,j11-j1+1,k11-k1+1);
            end
        end
        
        [~,fi] = min(cost);
        final_bone(j,k) = fi + ind_s;
        
        clear w;clear cost;clear tmp
    end
end


imag = zeros(x,y,z);
for j = 1:y
    for k = 1:z
        if ci(j,k) < 0.80
            imag(final_bone(j,k),j,k) = 1; % final bone surface
        end
    end
end
toc
save_untouch_nii(make_ana(imag,[.3 .3 .3], [0 0 0]),'Invivo2');
