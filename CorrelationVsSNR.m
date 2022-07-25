%% SNR-correlation comparisons
% note: this is between axons, not ROIs
cd('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\ACh\Preprocessed2P\SmallFiles\Combined')
snrfile = uigetfile('*.mat','MultiSelect','on');

corrs_all = [];
snrs_all = [];
for y = 1:length(snrfile)
    clearvars -except corrs_all snrs_all y snrfile
    load(char(snrfile(y)),'dFF_comb')
    
   for u = 1:size(dFF_comb,2)
       dFF_comb(:,u) = lowpass(dFF_comb(:,u),1,100);
   end
   
            %calculate signal to noise ratio
                resam_LP_F = dFF_comb';
                N = length(resam_LP_F); freq = 0:10/length(resam_LP_F):10/2;
            for i = 1:size(resam_LP_F,1)
                xdft(i,:) = fft(resam_LP_F(i,:)); 
                 
            end
                xdft = xdft(:,1:length(freq));
            for i = 1:size(resam_LP_F,1)
                psdx(i,:) = (1/(10*N)) * abs(xdft(i,:)).^2;
            end
                psdx(:,2:end-1) = 2*psdx(:,2:end-1);

                idx005 = find(freq>0.05);idx005 = idx005(1,1);
                idx05 = find(freq>0.5);idx05 = idx05(1,1);
                idx1 = find(freq>1);idx1 = idx1(1,1);
                idx3 = find(freq>3);idx3 = idx3(1,1);
            for i = 1:size(psdx,1)
                lowmax(i) = max(psdx(i,idx005:idx05));
                himax(i) = mean(psdx(i,idx1:idx3));
            end
                Ana_SNR = log(lowmax./himax);
                
                
      for t = 1:size(dFF_comb,2)
        for tt = 1:size(dFF_comb,2)
           %calculate correlation 
            temp = corrcoef(dFF_comb(:,t),dFF_comb(:,tt));
            corrs(t,tt) = temp(1,2);
            % make Ana_SNR variable the same size as corr
            SNRs(t,tt) = min([Ana_SNR(t),Ana_SNR(tt)]);
        end; end 
    
    
% remove one half of matrix and unfold
corrs = tril(corrs,-1);
SNRs = tril(SNRs,-1);

corrs(corrs==0) = NaN;
SNRs(SNRs==0)=NaN;

corrs = corrs(~isnan(corrs));
SNRs = SNRs(~isnan(SNRs));

corrs_all = [corrs_all corrs'];
snrs_all = [snrs_all SNRs'];
end

cd('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\ACh\Analysis\')
save('SNR-Corr_Exp_1LP.mat')

%% for same axon
cd('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\ACh\Preprocessed2P\SmallFiles\SameAxon')
snrfile = uigetfile('*.mat','MultiSelect','on');

corrs_all = [];
snrs_all = [];
for y = 1:length(snrfile)
    clearvars -except corrs_all snrs_all y snrfile
    load(char(snrfile(y)),'Fnorm')
    
   if size(Fnorm,1) == 2
            %calculate signal to noise ratio
   for u = 1:2
       Fnorm(u,:) = lowpass(Fnorm(u,:),1,100);
   end
   
                resam_LP_F = Fnorm;
                N = length(resam_LP_F); freq = 0:10/length(resam_LP_F):10/2;
            for i = 1:size(resam_LP_F,1)
                xdft(i,:) = fft(resam_LP_F(i,:));
            end
                xdft = xdft(:,1:length(freq));
            for i = 1:size(resam_LP_F,1)
                psdx(i,:) = (1/(10*N)) * abs(xdft(i,:)).^2;
            end
                psdx(:,2:end-1) = 2*psdx(:,2:end-1);

                idx005 = find(freq>0.05);idx005 = idx005(1,1);
                idx05 = find(freq>0.5);idx05 = idx05(1,1);
                idx1 = find(freq>1);idx1 = idx1(1,1);
                idx3 = find(freq>3);idx3 = idx3(1,1);
            for i = 1:size(psdx,1)
                lowmax(i) = max(psdx(i,idx005:idx05));
                himax(i) = mean(psdx(i,idx1:idx3));
            end
                Ana_SNR = log(lowmax./himax);
                
                
      for t = 1:size(Fnorm,1)
        for tt = 1:size(Fnorm,1) 
            % make Ana_SNR variable the same size as corr
            SNRs(t,tt) = min([Ana_SNR(t),Ana_SNR(tt)]);
        end; end 
    %calculate correlation 
            temp = corrcoef(Fnorm(1,:),Fnorm(2,:));
            corrs = temp(1,2);
    
    
% remove one half of matrix and unfold
SNRs = tril(SNRs,-1);
SNRs(SNRs==0)=NaN;
SNRs = SNRs(~isnan(SNRs));

corrs_all = [corrs_all corrs'];
snrs_all = [snrs_all SNRs'];
   end
end

corrs_all_same = corrs_all; 
snrs_all_same=snrs_all;

snrs_all_orig = snrs_all;
load('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\ACh\Analysis\SNR-Corr_Exp_1LP.mat')



figure; scatter(snrs_all,corrs_all)

hold on; scatter(snrs_all_same,corrs_all_same,'r')