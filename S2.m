function []= S2(Sessions);
prompt = {'Enter experiment type (ACh, NA, VCIN, or Widefield)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'Widefield'};
exptype = inputdlg(prompt,dlgtitle,dims,definput);
Dir1 = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\RawData\Spike2');
Dir2 = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Preprocessed2P\Spike2');
cd(Dir1)
Session=uigetfile('*.mat','MultiSelect','on');
Session = cellstr(Session);
for x = 1:length(Session)

  file = load(char(Session(x)));
    display('File loaded. Preprocessing should take about 2 minutes.')
    if strcmp(exptype, 'Widefield')
        V2P_fr_CL =  file.WF_Trig;
        
    else
V2P_fr_CL = file.V2P_fr_CL;
    end
Walk = file.Walk;
% if isfield(file,'Sound')
% Sound = file.Sound;
% end
% if isfield(file,'VNS')
% VNS = file.VNS;
% end

if isfield(file,'Lft_PPL')
    Pup=file.Lft_PPL;
end
if isfield(file,'Rt_PPL')
    Pup = file.Rt_PPL;
end

if isfield(file,'Lft_whsk')
    whisk=file.Lft_whsk;
end
if isfield(file,'Rt_whsk')
    whisk = file.Rt_whsk;
end

if isfield(file,'whisk')
    whisk=file.whisk;
end
if isfield(file,'Pup')
    Pup=file.Pup;
end  
if isfield(file,'Pupil')
    Pup=file.Pupil;
end 
if isfield(file,'Whisk')
    whisk = file.Whisk;
end
if isfield(file,'Whisking')
    whisk = file.Whisk;
end

% if isfield(file,'VNS')
% for y = 1:length(VNS.values)
%     if VNS.values(y) > 0.1
%         VNS.values(y) = 5;
%     else
%         VNS.values(y) = 0;
%     end
% end
% end

% if isfield(file,'Sound')
% for y = 1:length(Sound.values)
%     if Sound.values(y) > 0.5
%         Sound.values(y) = 5;
%     else Sound.values(y) = 0;
%     end
% end
% end

for y = 1:length(V2P_fr_CL.values)
    if V2P_fr_CL.values(y) > 0.5
        V2P_fr_CL.values(y) = 5;
    else
        V2P_fr_CL.values(y) = 0;
    end
end
Walk.values=(Walk.values-3)*10;
ts_encoder = timeseries(Walk.values,(linspace(0,(Walk.length*Walk.interval),Walk.length))');
% if isfield(file,'Sound')
% ts_sound = timeseries(Sound.values,(linspace(0,(Sound.length*Sound.interval),Sound.length))');
% else ts_sound = NaN;
% end
% if isfield(file,'VNS')
% ts_stim = timeseries(VNS.values,(linspace(0,(VNS.length*VNS.interval),VNS.length))');
% else ts_stim = NaN;
% end
ts_whisk = timeseries(whisk.values,(linspace(0,(whisk.length*whisk.interval),whisk.length))');
ts_2pFC = timeseries(V2P_fr_CL.values,(linspace(0,(V2P_fr_CL.length*V2P_fr_CL.interval),V2P_fr_CL.length))');
ts_pupil = timeseries(Pup.values,(linspace(0,(Pup.length*Pup.interval),Pup.length))');
figure ('Name','Pupil Raw Trace')
plot(ts_pupil.Data)
title('First choose maximum, then choose minimum pupil values.')
[~,y1] =(ginput(1));
[~,y2] = (ginput(1));
ts_pupil.Data(ts_pupil.Data < y2) = NaN;
ts_pupil.Data(ts_pupil.Data > y1) = NaN;
maxpup = max(ts_pupil.Data);
    ts_pupil.Data =(ts_pupil.Data/maxpup)*100;
    ts_pupil.Data=movmedian(ts_pupil.Data,2000);
%     figure('Name','Corrected Pupil Trace')
%     plot(ts_pupil.Data) 
 
%     histogram(ts_pupil.Data)

nanx = isnan(ts_pupil.Data);
t    = 1:numel(ts_pupil.Data);
ts_pupil.Data(nanx) = interp1(t(~nanx), ts_pupil.Data(~nanx), t(nanx));

% 
% if isfield(file,'VNS')
%  [vns_pks,vns_locs] = findpeaks(ts_stim.data,'MinPeakHeight',0.3, 'MinPeakDistance',1/VNS.interval*10);
% end
% figure ('Name','Stim locations');
% plot(ts_stim,ts_stim.time(vns_locs),vns_pks,'or')
% title('ts_stim peaks')
% if isfield(file,'Sound')
%  [sound_pks,sound_locs] = findpeaks(ts_sound.data,'MinPeakHeight',0.3, 'MinPeakDistance',1/Sound.interval*5);
%  figure ('Name','Sound locations');
%  plot(ts_sound,ts_sound.time(sound_locs),sound_pks,'or')
%  title('ts_sound peaks')
% end
FC_raw = ts_2pFC.data;
fc2p_locs_raw = find(FC_raw == 5);
fc2p_locs = find(FC_raw == 5);
for i = 1:length(fc2p_locs_raw)-1
    if fc2p_locs_raw(i+1)-fc2p_locs_raw(i)< 2
        fc2p_locs(i+1) =  NaN;
    end;end
fc2p_locs=fc2p_locs(~isnan(fc2p_locs));
fc2p_pks =  ones(length(fc2p_locs),1);
%figure; plot(FC_raw);
%hold on; plot(fc2p_locs,fc2p_pks,'o')
%[fc2p_pks,fc2p_locs] = findpeaks(ts_2pFC.data),'MinPeakHeight',0.3, 'MinPeakDistance',V2P_fr_CL.interval/2);
%figure ('Name','FrameClock locations');
%plot(ts_2pFC,ts_2pFC.time(fc2p_locs),fc2p_pks,'or')
%title('2P-Frame Clock Locations peaks')
FC_times=ts_2pFC.time(fc2p_locs);


% resample to 100hz
Resam_r = 100; Resam_p = 1/Resam_r;
times_resam = 0:Resam_p:max(ts_pupil.time);
tsresam_whisk = resample(ts_whisk,times_resam);
tsresam_walk = resample(ts_encoder, times_resam);
tsresam_2pFC = resample(ts_2pFC,times_resam);
tsresam_pupil = resample(ts_pupil,times_resam);
% if isfield(file,'Sound')
% tsresam_sound = resample(ts_sound,times_resam);
% else tsresam_sound = NaN;
% end
% if isfield(file,'VNS')
%     tsresam_stim = resample(ts_stim,times_resam);
% else tsresam_stim  =  NaN;
% end

%low pass filter to 10hz
resam_walk = tsresam_walk.data;
resam_pupil =  tsresam_pupil.data;
resam_whisk = tsresam_whisk.data;
% if isfield(file,'Sound')
% resam_sound = tsresam_sound.data;
% else resam_sound  = [];
% end
% if isfield(file,'VNS')
%     resam_stim = tsresam_stim.data;
% else resam_stim  = [];
% end
resam_pupil(isnan(resam_pupil)) = nanmean(resam_pupil);


resam_LP_walk = lowpass(resam_walk,10,100);
resam_LP_pupil = lowpass(resam_pupil,10,100);
resam_LP_whisk = lowpass(resam_whisk,10,100);
% if isfield(file,'Sound')
%     resam_LP_sound = lowpass(resam_sound,10,100);
% else resam_LP_sound = NaN;
%end
% if isfield(file,'VNS')
%     resam_LP_stim = lowpass(resam_stim,10,100);
% else resam_LP_stim = NaN;
% end


 cd(Dir2)
filename = strcat('S2_',char(Session(x)));
save(filename)
clearvars -except exptype Session Dir1 Dir2 x

disp('Done.')
close all
end
end

 