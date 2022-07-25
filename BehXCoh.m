%% BehCoh
clear
prompt = {'Enter experiment type (ACh or NA or VCIN)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'ACh'};
exptype = inputdlg(prompt,dlgtitle,dims,definput);

Dir1coh = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Preprocessed2P\SmallFiles\Combined');
Dir2coh = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Analysis\BehaviorCrossCorr');
cd(Dir1coh)
Sessions_coh=uigetfile('*.mat','Select the INPUT DATA FILE(s)','MultiSelect','on');
for q =  1:length(Sessions_coh)
    cd(Dir1coh)
    load(char(Sessions_coh(q)));
    
pupCoh=[];whiskCoh=[];walkCoh=[];
pupf=[];walkf=[];whiskf=[];


r = randperm(length(pupil));
pupil_shuf = pupil(r);
whisk_shuf = whisk(r);
walk_shuf = walk(r);


for i = 1:size(dFF_comb,2)


f = linspace(0.01,5);
 win=12000;
 ov=win*0.98;

 [pupCoh(:,i),pupf(:,i)] = mscohere(pupil,dFF_comb(:,i),hann(12000),12000*0.98,[],100);
  [whiskCoh(:,i),whiskf(:,i)] = mscohere(whisk,dFF_comb(:,i),hann(12000),12000*0.98,[],100);
   [walkCoh(:,i),walkf(:,i)] = mscohere(walk,dFF_comb(:,i),hann(12000),12000*0.98,[],100);
   
 [pupCoh_shuf(:,i),pupf_shuf(:,i)] = mscohere(pupil_shuf,dFF_comb(:,i),hann(12000),12000*0.98,[],100);
  [whiskCoh_shuf(:,i),whiskf_shuf(:,i)] = mscohere(whisk_shuf,dFF_comb(:,i),hann(12000),12000*0.98,[],100);
   [walkCoh_shuf(:,i),walkf_shuf(:,i)] = mscohere(walk_shuf,dFF_comb(:,i),hann(12000),12000*0.98,[],100);
end

% figure; stdshade(pupCoh',0.5,'m')
% hold on; stdshade(whiskCoh',0.5,'g')
% stdshade(walkCoh',0.5,'c')

%filename_coh = strcat('BehXCoh_',filename(7:end-10));
filename_coh = strcat('BehXCoh_shuf20220517_',filename(7:end-10));
clearvars -except pupCoh whiskCoh walkCoh pupf filename_coh q Sessions_coh Dir1coh Dir2coh pupCoh_shuf whiskCoh_shuf walkCoh_shuf
cd(Dir2coh)
save(filename_coh)
end


%% combining
prompt = {'Enter experiment type (ACh or NA or VCIN)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'ACh'};
exptype = inputdlg(prompt,dlgtitle,dims,definput);
% combining across all sessions
walk=[]; whisk=[];pupil=[];
walk_shuf=[]; whisk_shuf=[];pupil_shuf=[];
walk_lag=[];whisk_lag=[]; pupil_lag=[];
Dir = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Analysis\BehaviorCrossCorr');
cd(Dir)
Sessions=uigetfile('*.mat','Select the INPUT DATA FILE(s)','MultiSelect','on');
pupil=[];whisk=[];walk=[];f=[];
for q =  1:length(Sessions) 
    load(char(Sessions(q)));   
    pupil=[pupil pupCoh];
    whisk=[whisk whiskCoh];
    walk=[walk walkCoh];
        pupil_shuf=[pupil_shuf pupCoh_shuf];
    whisk_shuf=[whisk_shuf whiskCoh_shuf];
    walk_shuf=[walk_shuf walkCoh_shuf];
    f = [f pupf];
end
%
% figure; stdshade(pupil(1:34,:)',0.5,'m')
% hold on; stdshade(whisk(1:34,:)',0.5,'g')
% stdshade(walk(1:34,:)',0.5,'c')
% figure; plot(f(1:34,1),mean(whisk(1:34,:),2))

%below used to pliot like axon-axon
f_corrected = pupf(:,1)*(100/3.14);
figure; stdshade(pupil(5:411,:)',0.5,'m')
hold on; stdshade(whisk(5:411,:)',0.5,'g')
stdshade(walk(5:411,:)',0.5,'c')


figure; stdshade(pupil_shuf(5:411,:)',0.5,'m')
hold on; stdshade(whisk_shuf(5:411,:)',0.5,'g')
stdshade(walk_shuf(5:411,:)',0.5,'c')
