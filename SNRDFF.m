clear
disp('Preprocessing Suite2P output. This should take 3-5 minutes.')
cd('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project')
Metadata = readtable('AChNA_Meta.csv');

% Specify what session you'd like to work with and pull out relevant metadata
prompt = {'Enter experiment type (ACh or NA)', 'Enter Mouse ID', 'Enter date of recording','Enter session name'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'NA','3166','200721','3rois'};
Sess = inputdlg(prompt,dlgtitle,dims,definput);

% Pull out metadata to save with file
for o = 1:size(Metadata,1)
    temp(o,1) = strcmp(string(Metadata{o,2}),Sess(2,1)); 
    temp(o,2) = strcmp(string(Metadata{o,3}),Sess(3,1));
    temp(o,3) = strcmp(string(Metadata{o,4}),Sess(4,1));
    if sum(temp(o,:)) ==  3
        Meta_rows(o) = 1;
    else Meta_rows(o) =0;
    end
end
Meta_rows = find(Meta_rows);
for x = 1:length(Meta_rows)
meta_NumFrames(x) = table2array(Metadata(Meta_rows(x),8));
 meta_Session(x) = table2array(Metadata(Meta_rows(x),5));
 meta_NumZ(x) = table2array(Metadata(Meta_rows(x),6));
 meta_Fs(x) = table2array(Metadata(Meta_rows(x),7));
 meta_Plane(x) = table2array(Metadata(Meta_rows(x),9));
 meta_CenterX(x) = table2array(Metadata(Meta_rows(x),10));
 meta_CenterY(x) = table2array(Metadata(Meta_rows(x),11));
 meta_CenterZ(x) = table2array(Metadata(Meta_rows(x),12));
 meta_Ztop(x) = table2array(Metadata(Meta_rows(x),13));
  meta_AdjustX(x) = table2array(Metadata(Meta_rows(x),18));
 meta_AdjustY(x) = table2array(Metadata(Meta_rows(x),19));
 meta_ROIsize(x) = table2array(Metadata(Meta_rows(x),14));
end


%Dir1 = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(Sess(1,1)),'\RawData\Suite2P');
Dir1 = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(Sess(1,1)),'\RawData\Suite2P\SameAxon');
Recording = uigetdir(Dir1);

Rec_paths = genpath(Recording); 
Rec_C = strsplit(Rec_paths, ';');
Rec_dirs = []; 
for i = 1:length(Rec_C)  
    Rec_this_dir = string(Rec_C{i});   
    Rec_dirs = [Rec_dirs; Rec_this_dir];  
end

Rec_slash_counts = count(Rec_dirs,'\');
%Rec_files_idx = find(Rec_slash_counts == 10);
Rec_files_idx = find(Rec_slash_counts == 11);
Rec_planes = Rec_dirs(Rec_files_idx);


% start loop through planes
for u = 1:length(Rec_planes)
batch_this_Fall = load(strcat(Rec_planes(u),'\Fall.mat'));

%
F = [];
Fneu = [];

for j = 1:length(batch_this_Fall.iscell)
    if batch_this_Fall.iscell(j,1) ==1
F(j,:) = batch_this_Fall.F(j,:);
Fneu(j,:) = batch_this_Fall.Fneu(j,:);
    end
end
meta_NumFrames = size(F,2);
if size(F,1) > 0
F = nonzeros(F); Fneu = nonzeros(Fneu);
F = reshape(F,length(F)/meta_NumFrames,meta_NumFrames);
%Fneu = reshape(Fneu,length(Fneu)/meta_NumFrames,meta_NumFrames);


%pull in S2 file
cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(Sess(1,1)),'\Preprocessed2P\Spike2'));
load(strcat('S2_',char(meta_Session(1,1)),'.mat'))

% Now preprocess 2P data
 % Be sure you've already selected ROIs in Suite2P!
times_struct = FC_times(meta_NumZ:meta_NumZ:end);
if length(times_struct)>length(F)
times_struct = times_struct(1:size(F,2));
elseif length(F) > length(times_struct)
    F = F(:,1:length(times_struct));
end

 ts_F = timeseries(F(:,1:length(times_struct))',times_struct);
% upsample to 100hz
times_resam_2p = min(times_struct):Resam_p:max(times_struct);
tsresam_F = resample(ts_F,times_resam_2p);
%lowpass at 10 hz
resam_LP_F = lowpass(tsresam_F.data,10,Resam_r); resam_LP_F=resam_LP_F';
%calculate signal to noise ratio
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
if size(resam_LP_F,2) == length(resam_LP_pupil)
    resam_LP_F = resam_LP_F';
end


%exclude axons that don't meet SNR criteria
% for i = 1:size(resam_LP_F,1)
%     if Ana_SNR(i) < log(20)
%         resam_LP_F(i,:) = NaN;
%     end
% end
% resam_LP_F = resam_LP_F(~isnan(resam_LP_F));

 %meta_numKpt = size(resam_LP_F,1)/length(tsresam_F.data);
% 
% if meta_numKept >= 1
% resam_LP_F = reshape(resam_LP_F,meta_numKept,length(tsresam_F.data));
% else resam_LP_F = [];% end
meta_numKept = 1;

%calculate df/f
resam_LP_dff=[];
if meta_numKept >= 1
Fmed = median(resam_LP_F,2)';
for i = 1: size(resam_LP_F,1)
    resam_LP_dff(i,:) = (resam_LP_F(i,:) - Fmed(i))./Fmed(i);
end

% calculate Fnorm
for i = 1:size(resam_LP_F,1)
Fmin(i) = min(resam_LP_F(i,:));
end
for i = 1:size(resam_LP_F,1)
Fmax(i) = max(resam_LP_F(i,:));
end
resam_LP_Fnorm=[];
for i = 1:size(resam_LP_F,1)
    resam_LP_Fnorm(i,:) = (resam_LP_F(i,:) - Fmin(i))./(Fmax(i)-Fmin(i));
end
end

clearvars himax idx005 idx05 idx1 idx3 lowmax N psdx RL RU xdft freq 

% trim behavioral variables around 2p frame clock
On_2p = knnsearch(tsresam_pupil.time, times_resam_2p(1,1));
Off_2p = knnsearch(tsresam_pupil.time, times_resam_2p(1,end));
resam_LP_pupil = resam_LP_pupil(On_2p:Off_2p);
resam_LP_whisk = resam_LP_whisk(On_2p:Off_2p);
resam_LP_walk = resam_LP_walk(On_2p:Off_2p);
% if length(resam_LP_sound) > 2
% resam_LP_sound = resam_LP_sound(On_2p:Off_2p);
% end
% if length(resam_LP_stim) > 2
% resam_LP_stim = resam_LP_stim(On_2p:Off_2p);
% end
clearvars t*

%cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(Sess(1,1)),'\Preprocessed2P\Meso\'))
cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(Sess(1,1)),'\Preprocessed2P\Meso\SameAxon'))
filename = strcat('SNRDFF_',char(meta_Session(1)),'plane',num2str(meta_Plane(u)));
save(filename)
disp('Plane saved.')
end
end