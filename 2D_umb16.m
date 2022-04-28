%% ---------- Mohammad Arafat Hussain ---------------------%%
% comments and concerns to m.arafathussain@gmail.com
%
% ------------------------------------------------------------
clear all;
clc;
close all;

getd = @(p)path(p,path);
getd('C:\Users\arafat\Dropbox\ZPclustering');

[Im,header] = RPread_arafat('C:\Users\arafat\Dropbox\Elastography\06-10-2014-Small Parts\middle rib support.rf');
g = 21; %Best: 10, +3 
rf1_pre = Im(:,:,g);
rf2_pst = Im(:,:,g+4);

%% RF modification
rf1 = abs(hilbert(rf1_pre));
rf2 = abs(hilbert(rf2_pst));
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
btt = zeros(x,y);
ttb(1,:) = rf(1,:);
btt(x,:) = rf(x,:);

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
IRF = [-100 0];
alfa_DP = 0.2;
midA = 64;
alfa = 5;  
beta = 5;
gamma = 0.005;
T = 0.4;
f_0 = 10;
f_s = 40;
wDIff = 5;

str = am2d_mod_iv1(rf1,rf2,IRF,alfa_DP,midA,alfa,beta,gamma,T,f_0,f_s,wDIff);
%imagesc(str(:,1:end-20));axis off;axis square;colormap gray

%% strain processing
[a,b] = size(str);

for j = 1:b
    for i = 1:a
        if str(i,j) < 0
            str(i,j) = str(i,j)*(-1);
        end
    end
end
maxIm = max(str(:));
str = str/maxIm;
str = imadjust(str);
sn = str;

%% RF processing
rf = rf1;
rft = rf(43:size(rf,1)-43,11:size(rf,2)-10); % trimming to macth strain map
maxIm = max(rft(:));
rft = rft/maxIm;
rn = imadjust(rft);clear rft
% imagesc(rn(:,1:end-20));colormap gray;axis square;axis off

%% Automatic weight selection
n_init = 25;
svmin = 0;
svmax = 0.03;
str_init = svmin:(svmax-svmin)/(n_init-1):svmax;
crr_mat = zeros(1,n_init);

rf1_int = interpft(rf1_pre,6700); % BY A FACTOR OF 2.5
rf2_int = interpft(rf2_pst,6700); 

ROI1 = rf1_int(1:1500,51:80);clear rf1_int
ROI2 = rf2_int(1:1500,51:80);clear rf2_int
    
wsize = size(ROI1,1);
col = size(ROI2,2);

for k1 = 1:n_init
    seg_stch = (1-str_init(k1))*interp1q([1:wsize]',ROI2,(linspace(1,(1-str_init(k1))*(wsize),(wsize)))');
    for l = 1:col
        xx = ROI1(:,l);
        yy = seg_stch(:,l);
        [cc(:,l),lags] = xcorr(xx,yy,'coeff');
    end
    [C,I] = cosintp(mean(cc,2),lags);
    crr_mat(1,k1) = C;
end

ncc = max(crr_mat);

if ncc > 0.9
    w = 0.5;
elseif ncc < 0.5
    w = 0.1;
else
    w = ncc - 0.4;
end

%% Fused map
final = (1-w)*rn + w*sn;
% imagesc(final(:,1:end-20));colormap gray;axis square;axis off

%% Maximum points searching out
t = zeros(1,size(final,2));
tmp = max(final);
for i = 1:size(final,2)
     t(i) = max(find(tmp(i) == final(:,i)));
end

%% Deleting outlier
[~,l] = deleteoutliers(t);

if exist('l','var')
    t(l) = mean(t);
end

%% B-mode
env1 = envelope(rf1_pre);
option = 'mu';
fac = 255;
y1 = LogComp(env1,option,fac);
y11 = y1(43:size(rf,1)-43,11:size(rf,2)-30);

%% Kalman Filter
Data = ones(length(t),round(length(t)/4));
Data(:,1) = t';
[strain1 strain2] = Kal_LSQ(Data,13);
kb = strain1(:,1);
th = mean(abs(kb)) + 3*std(abs(kb));
m = 1;
h = 1;
flag = 0;
for i = 1:length(t)
    if h == 10;
        flag = 0;
    end
    if flag == 1
        h = h + 1;
        continue;
    end
    if abs(kb(i)) > th      
        loc(m) = i;
        m = m + 1;
        flag = 1;
    end
end

figure;imagesc(y11./max(y11(:)));
colormap gray;axis square;hold on;axis off;

%% Gaussian Mixture Regression 

if exist('loc','var')
   nbStates = length(loc) + 1;                % No. of segments
else
   nbStates = 1;
end

Data = [1:length(t);ones(1,length(t));t];
nbVar = size(Data,1); 
    
expData = zeros(1,length(t));
% figure; % Needed when ploting the GM
tr = 1;
for i = 1:nbStates
    if i == nbStates
        fi = length(t);
    else
        fi = loc(i) + 2;
    end
    
    % Training of GMM by EM algorithm, initialized by k-means clustering. 
    [NoGMM,~] = cluster(t(tr:fi),7);
    [Priors, Mu, Sigma] = EM_init_kmeans(Data(:,tr:fi), NoGMM);
    [Priors, Mu, Sigma] = EM(Data(:,tr:fi), Priors, Mu, Sigma);

    % GMR
    expData(tr:fi) = GMR(Priors,Mu,Sigma,tr:fi,1,nbVar);
    plot(tr:fi,expData(tr:fi),'r');hold on;
    tr = fi + 0;  
    
%     for n = 1:nbVar-1
%      subplot(1*(nbVar-1),1,1+(n-1)); hold on;
%      plotGMM(Mu([1,n+1],:), Sigma([1,n+1],[1,n+1],:), [0 .8 0], 1);
%      axis([min(Data(1,:)) max(Data(1,:)) min(Data(n+1,:))-0.01 max(Data(n+1,:))+0.01]);
%      xlabel('t','fontsize',16); ylabel(['x_' num2str(n)],'fontsize',16);
%     end
end