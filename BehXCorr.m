%% Cross correlations between axon data and behavioral data
% With time lag plots
% Q: How consistent is the timing of the relationship between behavior and axon activity? 
% FaceMap1_working has code for xcorr for facemap
% go from combined small file
function [file]= BehXCorr(Sessions_corr);
prompt = {'Enter experiment type (ACh or NA or VCIN)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'ACh'};
exptype = inputdlg(prompt,dlgtitle,dims,definput);

Dir1corr = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Preprocessed2P\SmallFiles\Combined');
Dir2corr = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Analysis\BehaviorCrossCorr');
cd(Dir1corr)
Sessions_corr=uigetfile('*.mat','Select the INPUT DATA FILE(s)','MultiSelect','on');
for q =  1:length(Sessions_corr)
    cd(Dir1corr)
    load(char(Sessions_corr(q)));

    pupil=lowpass(pupil,1,100);
    pupil_d = diff(pupil);
    whisk=lowpass(whisk,1,100);
    walk=lowpass(walk,1,100);
    for l = 1:size(Fnorm_comb,2)
    Fnorm_comb(:,l) = lowpass(dFF_comb(:,l),1,100);
end
    
Pupil_corr=[];Whisk_corr=[];Walk_corr=[];Pupil_d_corr=[];
for i = 1:size(Fnorm_comb,2)
this_corr = xcorr(pupil,Fnorm_comb(:,i),300,'coeff');
Pupil_corr(:,i) = this_corr;
this_corr = xcorr(pupil_d,Fnorm_comb(2:end,i),300,'coeff');
Pupil_d_corr(:,i) = this_corr;
this_corr = xcorr(whisk,Fnorm_comb(:,i),300,'coeff');
Whisk_corr(:,i) = this_corr;
this_corr = xcorr(walk,Fnorm_comb(:,i),300,'coeff');
Walk_corr(:,i) = this_corr;
if max(walk) < 0.05
    Walk_corr(:,i) = NaN;
end
end
for i = 1:size(Pupil_corr,2)
    Pupil_maxlag_corr(i) = find(Pupil_corr(:,i) == max(Pupil_corr(:,i)));
    Whisk_maxlag_corr(i) = find(Whisk_corr(:,i) == max(Whisk_corr(:,i)));
    Walk_maxlag_corr(i) = find(Walk_corr(:,i) == max(Walk_corr(:,i)));
end

mean_Pupil_corr = xcorr(pupil,mean(Fnorm_comb,2),300,'coeff');
mean_Pupil_d_corr = xcorr(pupil_d,mean(Fnorm_comb(2:end,:),2),300,'coeff');
mean_Whisk_corr = xcorr(whisk,mean(Fnorm_comb,2),300,'coeff');
mean_Walk_corr = xcorr(walk,mean(Fnorm_comb,2),300,'coeff');

% figure; subplot(1,3,1);stdshade(Pupil_corr'),title('Pupil')
% subplot(1,3,2);stdshade(Whisk_corr'),title('Whisk')
% subplot(1,3,3);stdshade(Walk_corr'),title('Walk')
 cd(Dir2corr)
% savefig(strcat('XCorrs_',filename(7:end-10),'.fig'))
% 
% figure; scatter(repmat(1,length(Pupil_maxlag_corr),1),Pupil_maxlag_corr)
% hold on;scatter(repmat(2,length(Pupil_maxlag_corr),1),Whisk_maxlag_corr)
% scatter(repmat(3,length(Pupil_maxlag_corr),1),Walk_maxlag_corr)
% xticks([1 2 3])
% xticklabels({'Pupil','Whisk','Walk'})
% savefig(strcat('ScatLags_',filename(7:end-10),'.fig'))

filename_corr = strcat('BehXCorr',filename(7:end-10));
clearvars -except *corr
clearvars this*
save(filename_corr)

end
end