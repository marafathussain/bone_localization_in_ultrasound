%%
function str = am2d_mod_iv1(rf1,rf2,IRF,alfa_DP,midA,alfa,beta,gamma,T,f_0,f_s,wDIff)


Im1 = rf1;
maxIm = max(Im1(:));
Im1 = Im1/maxIm;

Im2 = rf2;
Im2 = Im2/maxIm;

See_B_mode = 0;
if See_B_mode
    BMODE1 = log(abs(hilbert(Im1))+.01);
    figure,  imagesc(BMODE1(:,1:230));colormap(gray), colorbar; 
%     BMODE2 = log(abs(hilbert(Im2(40:end-40,10:end-10)))+.01);
%     figure, imagesc(BMODE2);colormap(gray), colorbar; axis square
end


% ---------------------------------------- %
% IRF = [-250 0];
IA = [-4 4]; %Maximum allowed disparity in lateral D
% alfa_DP = 0.15; % DP regularization weight
% ---------------------------------------- %
% ------------ 2D AM Paerametes ---------- %
% ---------------------------------------- %
% midA = 110;
% alfa = 25; %axial regularization
% beta = 5; %lateral regularization
% gamma = 0.005; %lateral regularization 
% T = .2; % threshold for IRLS
a_t = 1.00; % frequency dependent attenuation coefficient, in dB/cm/MHz
% f_0 = FileHeader.rfsd(1,1).TxFrequencyMhz/1e6; %ultrasound center freq. in MHz
% f_0=5;
% f_s = 20; % sampling freq. in MHz
xx = calc_att (a_t, f_0, f_s); % to compensate for attenuation

[D1 D2 DPdisparity] = AM2D(Im1, Im2, IRF, IA, midA, alfa_DP, alfa, beta, gamma, T, xx);
% the disp. of the first 40 and last 40 samples is not calculated in AM2D: 
% the disp. of the first 10 and last 10 A-lines is not calculated in AM2D: 
% figure, imagesc(D1), colorbar, title('axial displacement'), colormap(hot)

% ------------ Calculating Strain from Disp ---------- %
% ---------------------------------------------------- %
% wDIff = 43;
%size(D1)
[strain1 strain2] = Kal_LSQ(D1(41:end-41,11:end-10),wDIff);
%size(strain2)
strain1 = strain1((wDIff+1)/2:end-(wDIff-1)/2,:);
strain2 = strain2((wDIff+1)/2:end-(wDIff-1)/2,:);

%size(strain2)
startA = 1; endA = size(strain1,2);
startRF = 1; endRF = size(strain1,1); 
str = -strain2(startRF:endRF, startA:endA);
%Arafat
% [x1,y1] = size(str);
% 
% m = 0.2*max(str(:));
% 
% for rr = 1:x1
%     for hh = 1:y1
%         if str(rr,hh) >= m
%             str(rr,hh) = (-1)*str(rr,hh);
%         end
%     end
% end
% %
% dsp=D1;