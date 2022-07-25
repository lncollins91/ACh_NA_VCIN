%% Walking Onsets and Offsets
% From Combined Small File
function []= WalkOnOff(onoff_sessions);
prompt = {'Enter experiment type (ACh or NA)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'ACh'};
oType = inputdlg(prompt,dlgtitle,dims,definput);

cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(oType(1,1)),'\Preprocessed2P\SmallFiles\Combined'))
onoff_sessions = uigetfile('','MultiSelect','on'); 
for o = 1:length(onoff_sessions)
cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(oType(1,1)),'\Preprocessed2P\SmallFiles\Combined'))
clearvars -except onoff_sessions oType o
load(char(onoff_sessions(o)))

%
walk_bin = walk-median(walk);

for i = 1:length(walk_bin)  %criteria: 2.5 cm/s
    if walk_bin(i) > 0.025
        walk_bin(i) = 10;
    end;end
walk_flip=flip(walk_bin);

walk_onsets = find(walk_bin == 10);
walk_offsets = find(walk_flip == 10);
%
for i = 1:length(walk_onsets) % criteria: not within  2 s of beginning or end of recording
    if walk_onsets(i) < 301
        walk_onsets(i) = NaN;
    elseif walk_onsets(i) > length(walk)-301
        walk_onsets(i) = NaN;
    end;end;
for i = 1:length(walk_offsets)
    if walk_offsets(i) < 301
        walk_offsets(i) = NaN;
    elseif walk_offsets(i) > length(walk)-301
        walk_offsets(i) = NaN;
    end;end
walk_onsets = walk_onsets(~isnan(walk_onsets));
walk_offsets=walk_offsets(~isnan(walk_offsets));
% for i = 1:length(walk_onsets)
%     if mean(walk(walk_onsets(i):walk_onsets(i)+100)) < 0.01
%         walk_onsets(i) = NaN;
%     end
% end
% for i = 1:length(walk_offsets)
%     if mean(walk(walk_offsets(i) - 100:walk_offsets(i))) <0.01
%         walk_offsets(i) = NaN;
%     end
% end

walk_onsets = walk_onsets(~isnan(walk_onsets));
walk_offsets=walk_offsets(~isnan(walk_offsets));
%
temp_onsets = walk_onsets;
for i = 2:length(walk_onsets) % criteria: preceded by 1 s of stillness
    if temp_onsets(i)-temp_onsets(i-1) < 100
        walk_onsets(i) = NaN;
    end;end

temp_offsets = walk_offsets;
for i = 2:length(walk_offsets) % criteria: preceded by 1 s of stillness
    if temp_offsets(i)-temp_offsets(i-1) < 100
        walk_offsets(i) = NaN;
    end;end
walk_onsets = walk_onsets(~isnan(walk_onsets));
walk_offsets=walk_offsets(~isnan(walk_offsets));

%
walk_offsets = length(walk)-walk_offsets +1; % flip walk offsets back to normal time

for i = 1:length(walk_offsets)%% criteria: bout must last 1 s
    if walk_offsets(i) - walk_onsets(i) < 100
        walk_offsets(i) = NaN;
        walk_onsets(i) = NaN;
    end;end;
walk_onsets = walk_onsets(~isnan(walk_onsets));
walk_offsets=walk_offsets(~isnan(walk_offsets));
%
whisk_norm = round(((whisk-min(whisk))/(max(whisk)-min(whisk)))*100);
whisk_norm = whisk_norm-mode(whisk_norm);

for i = 1:length(walk_onsets)
    if walk_onsets(i) < 101
        walk_onsets(i) = NaN;
    end
    if walk_onsets(i) > length(pupil)-101
        walk_onsets(i) = NaN;
    end
end
walk_onsets = walk_onsets(~isnan(walk_onsets));
for i = 1:size(walk_onsets) %split traces for walk onset
    on_walking_whisk(:,i)= whisk_norm(walk_onsets(i) - 300:walk_onsets(i)+300);
    on_walking_pupil(:,i)= pupil(walk_onsets(i) - 300:walk_onsets(i)+300);
    on_walking_walk(:,i)= walk(walk_onsets(i) - 300:walk_onsets(i)+300);
    on_walking_dFF(:,:,i)= dFF_comb(walk_onsets(i) - 300:walk_onsets(i)+300,:);
    on_walking_Fnorm(:,:,i)= Fnorm_comb(walk_onsets(i) - 300:walk_onsets(i)+300,:);
end

for i = 1:size(walk_offsets) %split traces for walk offset
    off_walking_whisk(:,i)= whisk_norm(walk_offsets(i) - 300:walk_offsets(i)+300);
    off_walking_pupil(:,i)= pupil(walk_offsets(i) - 300:walk_offsets(i)+300);
    off_walking_walk(:,i)= walk(walk_offsets(i) - 300:walk_offsets(i)+300);
    off_walking_dFF(:,:,i)= dFF_comb(walk_offsets(i) - 300:walk_offsets(i)+300,:);
        off_walking_Fnorm(:,:,i)= Fnorm_comb(walk_offsets(i) - 300:walk_offsets(i)+300,:);
end
    
% 
% figure; subplot(2,2,1);stdshade(on_walking_walk');title('Walking Onset - Walk')
% subplot(2,2,2);stdshade(squeeze(mean(on_walking_dFF,3)'));title('Walking Onset - dFF')
% subplot(2,2,3);stdshade(off_walking_walk');title('Walking Offset - Walk')
% subplot(2,2,4);stdshade(squeeze(mean(off_walking_dFF,3)'));title('Walking Offset - dFF')


%
cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(oType(1,1)),'\Analysis\WalkOnOff'))
onoff_file = filename;
onTime_walk = walk_onsets;
offTime_walk=walk_offsets;
clearvars -except o*

save(strcat('Walk',onoff_file(6:end-10)))

end
%% to combine
prompt = {'Enter experiment type (ACh or NA)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'ACh'};
oType = inputdlg(prompt,dlgtitle,dims,definput);
cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(oType(1,1)),'\Analysis\WalkOnOff'))
comb_sessions = uigetfile('','MultiSelect','on'); 


Off_walking_dFF = [];
Off_walking_pupil = [];
Off_walking_whisk = [];
Off_walking_walk = [];
On_walking_dFF = [];
On_walking_pupil = [];
On_walking_whisk = [];
On_walking_walk = [];
On_walking_Fnorm=[];
Off_walking_Fnorm=[];


for i = 1:length(comb_sessions)
    load(char(comb_sessions(i)))
    
if exist('off_walking_dFF')
Off_walking_dFF = [Off_walking_dFF squeeze(mean(off_walking_dFF,3))];
Off_walking_Fnorm = [Off_walking_Fnorm squeeze(mean(off_walking_Fnorm,3))];
Off_walking_pupil = [Off_walking_pupil off_walking_pupil];
Off_walking_whisk = [Off_walking_whisk off_walking_whisk];
Off_walking_walk = [Off_walking_walk off_walking_walk];
end
if exist('on_walking_dFF')
On_walking_dFF = [On_walking_dFF squeeze(mean(on_walking_dFF,3))];
On_walking_Fnorm = [On_walking_Fnorm squeeze(mean(on_walking_Fnorm,3))];
On_walking_pupil = [On_walking_pupil on_walking_pupil];
On_walking_whisk = [On_walking_whisk on_walking_whisk];
On_walking_walk = [On_walking_walk on_walking_walk]; 
end

COUNT = size(Off_walking_dFF,2);

end


figure; subplot(2,2,1);stdshade(On_walking_walk');title('Walking Onset - Walk');ylim([0 0.12])
subplot(2,2,2);stdshade(Off_walking_walk');title('Walking Offset - Walk');ylim([0 0.12])
subplot(2,2,3);stdshade(squeeze(mean(On_walking_dFF,3)'));title('Walking Onset - dFF');ylim([-0.1 3])
subplot(2,2,4);stdshade(squeeze(mean(Off_walking_dFF,3)'));title('Walking Offset - dFF');ylim([-0.1,3])

clearvars -except O*

save('Combined_WalkOnOff')
%save('Combined_WalkOnOff_bleb')
end


