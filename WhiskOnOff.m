%% Whisking Onsets and Offsets
% From Combined Small File
function []= WhiskOnOff(onoff_sessions);
prompt = {'Enter experiment type (ACh or NA or Widefield)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'ACh'};
oType = inputdlg(prompt,dlgtitle,dims,definput);

cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(oType(1,1)),'\Preprocessed2P\SmallFiles\SameAxon'))
onoff_sessions = uigetfile('','MultiSelect','on'); 
for o = 1:length(onoff_sessions)
cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(oType(1,1)),'\Preprocessed2P\SmallFiles\SameAxon'))
clearvars -except onoff_sessions oType o
load(char(onoff_sessions(o)))

%% 
whisk_norm = round(((whisk-min(whisk))/(max(whisk)-min(whisk)))*100);
whisk_norm = whisk_norm-mode(whisk_norm);
whisk_flip = flip(whisk_norm);


% locs = findpeaks(whisk_norm);
% whisk_onsets=cell2mat(struct2cell(locs));
[pks,locs] = findpeaks(whisk_norm);
whisk_onsets=locs;
% locs =  findpeaks(whisk_flip);
[pks, locs] = findpeaks(whisk_flip);
% whisk_offsets=cell2mat(struct2cell(locs));
whisk_offsets = locs;

for i = 1:length(whisk_onsets) %clip off onsets within 5 s of beginning or end of recording
    if whisk_onsets(i) < 500
        whisk_onsets(i) = NaN;
    elseif whisk_onsets(i) > length(whisk)-500
        whisk_onsets(i) = NaN;
    end;end
whisk_onsets = whisk_onsets(~isnan(whisk_onsets));
for i = 1:length(whisk_offsets) %clip off onsets within 5 s of beginning or end of recording
    if whisk_offsets(i) < 500
        whisk_offsets(i) = NaN;
    elseif whisk_offsets(i) > length(whisk)-500
        whisk_offsets(i) = NaN;
    end;end
whisk_offsets = whisk_offsets(~isnan(whisk_offsets));

%%
for i=  1:length(whisk_onsets)%remove whisking periods with < 20% normalized whisk
    if whisk_norm(whisk_onsets(i)) < 20
        whisk_onsets(i) = NaN;
    end;end;
whisk_onsets = whisk_onsets(~isnan(whisk_onsets));
for i=  1:length(whisk_offsets)%remove whisking periods with < 20% normalized whisk
    if whisk_flip(whisk_offsets(i)) < 20
        whisk_offsets(i) = NaN;
    end;end;
whisk_offsets = whisk_offsets(~isnan(whisk_offsets));

%%
temp_onsets = whisk_onsets;
temp_offsets = whisk_offsets;
for i =  2:length(whisk_onsets) % require 1 s silence before onset
    if temp_onsets(i) - temp_onsets(i-1) < 100
    whisk_onsets(i) = NaN;
    end;end
whisk_onsets = whisk_onsets(~isnan(whisk_onsets));
for i =  2:length(whisk_offsets) % require 1 s silence before onset
    if temp_offsets(i) - temp_offsets(i-1) < 100
    whisk_offsets(i) = NaN;
    end;end
whisk_offsets = whisk_offsets(~isnan(whisk_offsets));
%%

twitch_onsets = whisk_onsets;
twitch_offsets = whisk_offsets;
for i = 1:length(whisk_onsets) % require 1 s of movement for whisking
    if mean(whisk_norm(whisk_onsets(i):whisk_onsets(i)+100)) < 10
        whisk_onsets(i) = NaN;
    else twitch_onsets(i) = NaN;
    end;end;
whisk_onsets = whisk_onsets(~isnan(whisk_onsets));
twitch_onsets = twitch_onsets(~isnan(twitch_onsets));
for i = 1:length(whisk_offsets) % require 1 s of movement for whisking
    if mean(whisk_flip(whisk_offsets(i):whisk_offsets(i)+100)) < 10
        whisk_offsets(i) = NaN;
    else twitch_offsets(i) = NaN;
    end;end;
whisk_offsets = whisk_offsets(~isnan(whisk_offsets));
twitch_offsets = twitch_offsets(~isnan(twitch_offsets));
%%
% figure; scatter(whisk_onsets,zeros(length(whisk_onsets),1));hold on;
% scatter(twitch_onsets,zeros(length(twitch_onsets),1));
% plot(whisk_norm)

for i = 1:size(whisk_onsets) %split traces for whisk onset
    on_whisking_whisk(:,i)= whisk_norm(whisk_onsets(i) - 300:whisk_onsets(i)+300);
    on_whisking_pupil(:,i)= pupil(whisk_onsets(i) - 300:whisk_onsets(i)+300);
    on_whisking_walk(:,i)= walk(whisk_onsets(i) - 300:whisk_onsets(i)+300);
    on_whisking_dFF(:,:,i)= dFF_comb(whisk_onsets(i) - 300:whisk_onsets(i)+300,:);
    on_whisking_Fnorm(:,:,i) = Fnorm_comb(whisk_onsets(i)- 300:whisk_onsets(i)+300,:);
end

for i = 1:size(twitch_onsets) %split traces for twitch onset
    on_twitch_whisk(:,i)= whisk_norm(twitch_onsets(i) - 300:twitch_onsets(i)+300);
    on_twitch_pupil(:,i)= pupil(twitch_onsets(i) - 300:twitch_onsets(i)+300);
    on_twitch_walk(:,i)= walk(twitch_onsets(i) - 300:twitch_onsets(i)+300);
    on_twitch_dFF(:,:,i)= dFF_comb(twitch_onsets(i) - 300:twitch_onsets(i)+300,:);
    on_twitch_Fnorm(:,:,i)= Fnorm_comb(twitch_onsets(i) - 300:twitch_onsets(i)+300,:);
end

pupil_flip=flip(pupil);walk_flip=flip(walk);
for i = 1:size(dFF_comb,2)
    dFF_flip(:,i) = flip(dFF_comb(:,i));
    Fnorm_flip(:,i) = flip(Fnorm_comb(:,i));
end
for i = 1:size(whisk_offsets) %split traces for whisk offset
    off_whisking_whisk(:,i)= whisk_flip(whisk_offsets(i) - 300:whisk_offsets(i)+300);
    off_whisking_pupil(:,i)= pupil_flip(whisk_offsets(i) - 300:whisk_offsets(i)+300);
    off_whisking_walk(:,i)= walk_flip(whisk_offsets(i) - 300:whisk_offsets(i)+300);
    off_whisking_dFF(:,:,i)= dFF_flip(whisk_offsets(i) - 300:whisk_offsets(i)+300,:);
    off_whisking_Fnorm(:,:,i)= Fnorm_flip(whisk_offsets(i) - 300:whisk_offsets(i)+300,:);
end

for i = 1:size(twitch_offsets) %split traces for twitch offset
    off_twitch_whisk(:,i)= whisk_flip(twitch_offsets(i) - 300:twitch_offsets(i)+300);
    off_twitch_pupil(:,i)= pupil_flip(twitch_offsets(i) - 300:twitch_offsets(i)+300);
    off_twitch_walk(:,i)= walk_flip(twitch_offsets(i) - 300:twitch_offsets(i)+300);
    off_twitch_dFF(:,:,i)= dFF_flip(twitch_offsets(i) - 300:twitch_offsets(i)+300,:);
    off_twitch_Fnorm(:,:,i)= Fnorm_flip(twitch_offsets(i) - 300:twitch_offsets(i)+300,:);
end

%%
% flip offsets back to normal time

off_whisking_whisk = flip(off_whisking_whisk);
off_whisking_pupil = flip(off_whisking_pupil);
off_whisking_walk = flip(off_whisking_walk);
off_twitch_whisk = flip(off_twitch_whisk);
off_twitch_pupil = flip(off_twitch_pupil);
off_twitch_walk = flip(off_twitch_walk);
for i = 1:size(off_whisking_dFF,2)
    off_whisking_dFF(:,i,:) = flip(off_whisking_dFF(:,i,:));
    off_whisking_Fnorm(:,i,:) = flip(off_whisking_Fnorm(:,i,:));
end
for i = 1:size(off_twitch_dFF,2)
    off_twitch_dFF(:,i,:) = flip(off_twitch_dFF(:,i,:));
    off_twitch_Fnorm(:,i,:) = flip(off_twitch_Fnorm(:,i,:));
end
    

figure; subplot(2,4,1);stdshade(on_whisking_whisk');title('Whisking Onset - Whisk')
subplot(2,4,3);stdshade(squeeze(mean(on_whisking_dFF,3)'));title('Whisking Onset - dFF')
subplot(2,4,2);stdshade(on_twitch_whisk');title('Twitch Onset - Whisk')
subplot(2,4,4);stdshade(squeeze(mean(on_twitch_dFF,3)'));title('Twitch Onset - dFF')
subplot(2,4,5);stdshade(off_whisking_whisk');title('Whisking Offset - Whisk')
subplot(2,4,7);stdshade(squeeze(mean(off_whisking_dFF,3)'));title('Whisking Offset - dFF')
subplot(2,4,6);stdshade(off_twitch_whisk');title('Twitch Offset - Whisk')
subplot(2,4,8);stdshade(squeeze(mean(off_twitch_dFF,3)'));title('Twitch Offset - dFF')


%%
cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(oType(1,1)),'\Analysis\WhiskOnOff'))
onoff_file = filename;
onTime_whisk = whisk_onsets;
offTime_whisk=whisk_offsets;
onTime_twitch=twitch_onsets;
offTime_twitch=twitch_offsets;
clearvars -except o*

change_whisking_whisk = max(on_whisking_whisk(100:200,:))-mean(on_whisking_whisk(50:100,:));
change_whisking_walk = max(on_whisking_walk(100:200,:))-mean(on_whisking_walk(50:100,:));
change_whisking_pupil = max(on_whisking_pupil(100:200,:))-mean(on_whisking_pupil(50:100,:));
for i = 1:size(on_whisking_dFF,2)
change_whisking_dFF(i,:) = max(on_whisking_dFF(100:200,i,:))-mean(on_whisking_dFF(50:100,i,:));
change_whisking_Fnorm(i,:) = max(on_whisking_Fnorm(100:200,i,:))-mean(on_whisking_Fnorm(50:100,i,:));
end
change_twitch_whisk = max(on_twitch_whisk(100:200,:))-mean(on_twitch_whisk(50:100,:));
change_twitch_walk = max(on_twitch_walk(100:200,:))-mean(on_twitch_walk(50:100,:));
change_twitch_pupil = max(on_twitch_pupil(100:200,:))-mean(on_twitch_pupil(50:100,:));
for i = 1:size(on_twitch_dFF,2)
change_twitch_dFF(i,:) = max(on_twitch_dFF(100:200,i,:))-mean(on_twitch_dFF(50:100,i,:));
change_twitch_Fnorm(i,:) = max(on_twitch_Fnorm(100:200,i,:))-mean(on_twitch_Fnorm(50:100,i,:));
end

save(strcat('Whisk',onoff_file(6:end-10)))

end
%% START HERE TO JUST MAKE COMBINED FIGS

prompt = {'Enter experiment type (ACh or NA or Widefield)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'ACh'};
oType = inputdlg(prompt,dlgtitle,dims,definput);

cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(oType(1,1)),'\Analysis\WhiskOnOff'))
comb_sessions = uigetfile('','MultiSelect','on'); 

Off_twitch_dFF = [];
Off_twitch_pupil = [];
Off_twitch_whisk = [];
Off_twitch_walk = [];
On_twitch_dFF = [];
On_twitch_pupil = [];
On_twitch_whisk = [];
On_twitch_walk = [];
Off_whisking_dFF = [];
Off_whisking_pupil = [];
Off_whisking_whisk = [];
Off_whisking_walk = [];
On_whisking_dFF = [];
On_whisking_pupil = [];
On_whisking_whisk = [];
On_whisking_walk = [];
Off_twitch_Fnorm=[];
On_twitch_Fnorm=[];
On_whisking_Fnorm=[];
Off_whisking_Fnorm=[];


for i = 1:length(comb_sessions)
    load(char(comb_sessions(i)))
    
Off_twitch_dFF = [Off_twitch_dFF squeeze(mean(off_twitch_dFF,3))];
Off_twitch_pupil = [Off_twitch_pupil off_twitch_pupil];
Off_twitch_whisk = [Off_twitch_whisk off_twitch_whisk];
Off_twitch_walk = [Off_twitch_walk off_twitch_walk];
On_twitch_dFF = [On_twitch_dFF squeeze(mean(on_twitch_dFF,3))];
On_twitch_pupil = [On_twitch_pupil on_twitch_pupil];
On_twitch_whisk = [On_twitch_whisk on_twitch_whisk];
On_twitch_walk = [On_twitch_walk on_twitch_walk];
Off_whisking_dFF = [Off_whisking_dFF squeeze(mean(off_whisking_dFF,3))];
Off_whisking_pupil = [Off_whisking_pupil off_whisking_pupil];
Off_whisking_whisk = [Off_whisking_whisk off_whisking_whisk];
Off_whisking_walk = [Off_whisking_walk off_whisking_walk];
On_whisking_dFF = [On_whisking_dFF squeeze(mean(on_whisking_dFF,3))];
On_whisking_pupil = [On_whisking_pupil on_whisking_pupil];
On_whisking_whisk = [On_whisking_whisk on_whisking_whisk];
On_whisking_walk = [On_whisking_walk on_whisking_walk]; 

On_whisking_Fnorm = [On_whisking_Fnorm squeeze(mean(on_whisking_Fnorm,3))];
Off_whisking_Fnorm = [Off_whisking_Fnorm squeeze(mean(off_whisking_Fnorm,3))];
On_twitch_Fnorm = [On_twitch_Fnorm squeeze(mean(on_twitch_Fnorm,3))];
Off_twitch_Fnorm = [Off_twitch_Fnorm squeeze(mean(off_twitch_Fnorm,3))];
end

% figure; subplot(2,4,1);stdshade(On_whisking_whisk');title('Whisking Onset - Whisk');ylim([0 40])
% subplot(2,4,2);stdshade(Off_whisking_whisk');title('Whisking Offset - Whisk');ylim([0 40])
% subplot(2,4,3);stdshade(On_twitch_whisk');title('Twitch Onset - Whisk');ylim([0 40])
% subplot(2,4,4);stdshade(Off_twitch_whisk');title('Twitch Offset - Whisk');ylim([0 40])
% subplot(2,4,5);stdshade(squeeze(mean(On_whisking_dFF,3)'));ylim([-0.1 0.3]);title('Whisking Onset - dFF');ylim([-0.1 1])
% subplot(2,4,6);stdshade(squeeze(mean(Off_whisking_dFF,3)'));ylim([-0.1 0.3]);title('Whisking Offset - dFF');ylim([-0.1 1])
% subplot(2,4,7);stdshade(squeeze(mean(On_twitch_dFF,3)'));ylim([-0.1 0.3]);title('Twitch Onset - dFF');ylim([-0.1 1])
% subplot(2,4,8);stdshade(squeeze(mean(Off_twitch_dFF,3)'));ylim([-0.1 0.3]);title('Twitch Offset - dFF');ylim([-0.1 1])



save('Combined_WhiskOnOff')

% figure;
% for i = 1:length(comb_sessions)
%     load(char(comb_sessions(i)))
%    scatter(change_whisking_whisk,mean(change_whisking_dFF),[],'red')
%    hold on; scatter(change_twitch_whisk,mean(change_twitch_dFF),[],'black');
% xlabel('Whisker Pad MEI'); ylabel('dF/F')
% 
figure; subplot(2,4,1);stdshade(On_whisking_whisk');title('Whisking Onset - Whisk')
subplot(2,4,2);stdshade(Off_whisking_whisk');title('Whisking Offset - Whisk')
subplot(2,4,3);stdshade(On_twitch_whisk');title('Twitch Onset - Whisk')
subplot(2,4,4);stdshade(Off_twitch_whisk');title('Twitch Offset - Whisk')
subplot(2,4,5);stdshade(On_whisking_dFF');title('Whisking Onset - dFF')
subplot(2,4,6);stdshade(Off_whisking_dFF');title('Whisking Offset - dFF')
subplot(2,4,7);stdshade(On_twitch_dFF');title('Twitch Onset - dFF')
subplot(2,4,8);stdshade(Off_twitch_dFF');title('Twitch Offset - dFF')

% confidence interval 

% figure; 
% for i = 1:size(On_whisking_dFF,2)
%     hold on;
%     plot(On_whisking_dFF(:,i)+i)
% end


%end
end


