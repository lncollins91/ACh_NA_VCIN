%% Same axon coherence to see if coherence is high at low frequencies 
% to explain the low correlation with a 1s moving window
% 4/5/22 this has been repurposed to be Fig4 D and E.

%% Correlations: overall and Whisk and Walk onset/offset 
% make a SameAxon folder in Analysis\WindowedCorrelation
% below is copied from WinCorr:
% this is used for Fig 4
clear
prompt = {'Enter experiment type (ACh or NA or VCIN)'};
dlgtitle = 'Input';
dims = [1 35];

definput = {'ACh'};
exptype = inputdlg(prompt,dlgtitle,dims,definput);

Dir1coh = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Preprocessed2P\SmallFiles\SameAxon');
Dir2coh = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Analysis\WindowedCorrelation\SameAxon');


overallcoh = [];


cd(Dir1coh)
SessCoh=uigetfile('*.mat','MultiSelect','on');
for x =1:length(SessCoh)
    clearvars -except SessCoh Dir1coh Dir2coh exptype x  overallcoh 
    load(char(SessCoh(x)),'Fnorm','filename');

 f = linspace(0.01,5);
 win=12000
 ov=win*0.98;
[allcoh_2mWin,fcoh_2mwin] = mscohere(Fnorm(1,:),Fnorm(2,:),win,ov,f,100);
 win=6000;
 ov=win*0.98;
[allcoh_1mWin,fcoh_1mwin] = mscohere(Fnorm(1,:),Fnorm(2,:),win,ov,f,100);
 win=3000;
 ov=win*0.98;
[allcoh_30sWin,fcoh_30swin] = mscohere(Fnorm(1,:),Fnorm(2,:),win,ov,f,100);
 win=1000;
 ov=win*0.98;
[allcoh_10sWin,fcoh_10swin] = mscohere(Fnorm(1,:),Fnorm(2,:),win,ov,f,100);
 win=500;
 ov=win*0.98;
[allcoh_5sWin,fcoh_5swin] = mscohere(Fnorm(1,:),Fnorm(2,:),win,ov,f,100);
%figure; plot(f,allcoh)

cd(Dir2coh)

save(strcat('sameCoh_MultWin_',filename(7:end-4)))
    %overallcoh=[overallcoh allcoh];
end

%% combining sameCoh 
Dir2coh = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\NA\Analysis\WindowedCorrelation\SameAxon');

overallcoh_2mWin = [];
overallcoh_1mWin = [];
overallcoh_30sWin = [];
overallcoh_10sWin = [];
overallcoh_5sWin = [];

cd(Dir2coh)
Sss=uigetfile('*.mat','MultiSelect','on');
for xx =1:length(Sss)
    clearvars -except Sss  Dir2coh  xx  overallcoh_2mWin overallcoh_1mWin overallcoh_30sWin overallcoh_10sWin overallcoh_5sWin  
    load(char(Sss(xx)));
overallcoh_2mWin = [overallcoh_2mWin allcoh_2mWin'];
overallcoh_1mWin = [overallcoh_1mWin allcoh_1mWin'];
overallcoh_30sWin = [overallcoh_30sWin allcoh_30sWin'];
overallcoh_10sWin = [overallcoh_10sWin allcoh_10sWin'];
overallcoh_5sWin = [overallcoh_5sWin allcoh_5sWin'];

end
save('SameCoh_MultWin_all')
%%

% figure; stdshade(overallcoh_2m');title('2mwin')
% figure; stdshade(overallcoh_1m');title('1mwin')
% figure; stdshade(overallcoh_30s');title('30swin')
% figure; stdshade(overallcoh_10s');title('10swin')
% figure; stdshade(overallcoh_5s');title('5swin')


figure; stdshade(overallcoh_2mWin',0.5,'k');
hold on; stdshade(overallcoh_1mWin',0.5,'b');
stdshade(overallcoh_30sWin',0.5,'g');
stdshade(overallcoh_10sWin',0.5,'y');
stdshade(overallcoh_5sWin',0.5,'m');

%% Now repeat exactly for differnt axon data
clear
prompt = {'Enter experiment type (ACh or NA or VCIN)'};
dlgtitle = 'Input';
dims = [1 35];

definput = {'ACh'};
exptype = inputdlg(prompt,dlgtitle,dims,definput);

Dir1coh = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Preprocessed2P\SmallFiles\Combined');
Dir2coh = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Analysis\WindowedCorrelation');


overallcoh = [];


cd(Dir1coh)
SessCoh=uigetfile('*.mat','MultiSelect','on');
for x =1:length(SessCoh)
    clearvars -except SessCoh Dir1coh Dir2coh exptype x  overallcoh
    load(char(SessCoh(x)),'Fnorm_comb','filename');

    Fnorm=Fnorm_comb';
 f = linspace(0.01,5);
 win=12000;
 ov=win*0.98;
[allcoh_2mWin,fcoh_2mwin] = mscohere(Fnorm(1,:),Fnorm(2,:),win,ov,f,100);
 win=6000;
 ov=win*0.98;
[allcoh_1mWin,fcoh_1mwin] = mscohere(Fnorm(1,:),Fnorm(2,:),win,ov,f,100);
 win=3000;
 ov=win*0.98;
[allcoh_30sWin,fcoh_30swin] = mscohere(Fnorm(1,:),Fnorm(2,:),win,ov,f,100);
 win=1000;
 ov=win*0.98;
[allcoh_10sWin,fcoh_10swin] = mscohere(Fnorm(1,:),Fnorm(2,:),win,ov,f,100);
 win=500;
 ov=win*0.98;
[allcoh_5sWin,fcoh_5swin] = mscohere(Fnorm(1,:),Fnorm(2,:),win,ov,f,100);
%figure; plot(f,allcoh)

cd(Dir2coh)

save(strcat('difCoh_MultWin_',filename(7:end-4)))
    %overallcoh=[overallcoh allcoh];
end

% combining difCoh 
Dir2coh = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\NA\Analysis\WindowedCorrelation');

overallcoh_2mWin = [];
overallcoh_1mWin = [];
overallcoh_30sWin = [];
overallcoh_10sWin = [];
overallcoh_5sWin = [];

cd(Dir2coh)
Sss=uigetfile('*.mat','MultiSelect','on');
for xx =1:length(Sss)
    clearvars -except Sss  Dir2coh  xx  overallcoh_2mWin overallcoh_1mWin overallcoh_30sWin overallcoh_10sWin overallcoh_5sWin  
    load(char(Sss(xx)));
overallcoh_2mWin = [overallcoh_2mWin allcoh_2mWin'];
overallcoh_1mWin = [overallcoh_1mWin allcoh_1mWin'];
overallcoh_30sWin = [overallcoh_30sWin allcoh_30sWin'];
overallcoh_10sWin = [overallcoh_10sWin allcoh_10sWin'];
overallcoh_5sWin = [overallcoh_5sWin allcoh_5sWin'];

end
save('difCoh_MultWin_all')%%

% figure; stdshade(overallcoh_2m');title('2mwin')
% figure; stdshade(overallcoh_1m');title('1mwin')
% figure; stdshade(overallcoh_30s');title('30swin')
% figure; stdshade(overallcoh_10s');title('10swin')
% figure; stdshade(overallcoh_5s');title('5swin')


figure; stdshade(overallcoh_2mWin',0.5,'k');set(gca, 'XScale', 'log')
hold on; stdshade(overallcoh_1mWin',0.5,'b');
stdshade(overallcoh_30sWin',0.5,'g');
stdshade(overallcoh_10sWin',0.5,'y');
stdshade(overallcoh_5sWin',0.5,'m');

%% Now repeat exactly for differnt axon data axon-axon, not ROI-ROI
clear
prompt = {'Enter experiment type (ACh or NA or VCIN)'};
dlgtitle = 'Input';
dims = [1 35];

definput = {'ACh'};
exptype = inputdlg(prompt,dlgtitle,dims,definput);

Dir1coh = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Preprocessed2P\SmallFiles\Combined');
Dir2coh = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Analysis\WindowedCorrelation');

overallcoh_2mWin = [];


cd(Dir1coh)
SessCoh=uigetfile('*.mat','MultiSelect','on');
for x =1:length(SessCoh)
    cd(Dir1coh)
    clearvars -except SessCoh Dir1coh Dir2coh exptype x overallcoh_2mWin
    load(char(SessCoh(x)),'Fnorm_comb','filename');

%     if size(Fnorm_comb,2) > 14
%         Fnorm_comb = Fnorm_comb(:,1:15);
%     end
    

Fnorm = lowpass(Fnorm_comb,10,100);
Fnorm = Fnorm';

 f = linspace(0.01,5);
 win=18000;
 ov=win*0.98;
%if size(Fnorm,1) > 4
    allcoh_2mWin=[];temp=[];
 for t = 1:size(Fnorm,1)
     for tt = 1:size(Fnorm,1)
%[temp,fcoh] = mscohere(Fnorm(t,:),Fnorm(tt,:),win,ov,f,100);
%below is updated version from 4/27/22
% nfft = 12000; %2m window used to calculate fft
%  noverlap =nfft*0.98; %98 percent overlap
%  [coh,fcoh] = mscohere(Fnorm(t,:),Fnorm(tt,:),hann(nfft),noverlap);
% allcoh_2mWin(t,tt,:) = coh;
%below to make final plots(as of 5/4/22)
 [coh,fcoh] = mscohere(Fnorm(t,:),Fnorm(tt,:),hann(12000),12000*0.98,[],100);
 allcoh_2mWin(t,tt,:)=coh;
 
    % end
 end
end
for u = 1:100;
allcoh_2mWin_tril(:,:,u) = tril(squeeze(allcoh_2mWin(:,:,u)),-1);
end
allcoh_2mWin_tril(allcoh_2mWin_tril==0) = NaN;
allcoh_2mWin_tril=nanmean(allcoh_2mWin_tril,1);
allcoh_2mWin_tril=nanmean(allcoh_2mWin_tril,1);
allcoh_2mWin_tril=squeeze(allcoh_2mWin_tril);


overallcoh_2mWin = [overallcoh_2mWin allcoh_2mWin_tril'];
cd(Dir2coh)

save(strcat('difCoh_axonaxon_final20220504_10hzLP',filename(7:end-4)))
end

%save('difCoh_axonaxon_2mWin_final20220504')


%% combine axon-axon
clear
prompt = {'Enter experiment type (ACh or NA or VCIN)'};
dlgtitle = 'Input';
dims = [1 35];

definput = {'ACh'};
exptype = inputdlg(prompt,dlgtitle,dims,definput);

cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Analysis\WindowedCorrelation'));
figfiles = uigetfile('*.mat','MultiSelect','on');
all=[];
all_sep=[];

for f = 1:length(figfiles)
    load(char(figfiles(f)),'allcoh_2mWin_tril','fcoh','allcoh_2mWin')
    atril=[];
    for u = 1:821;
atril(:,:,u) = tril(squeeze(allcoh_2mWin(:,:,u)),-1);
    end
atril(atril==0) = NaN;
atril=nanmean(atril,1);
atril=nanmean(atril,1);
atril=squeeze(atril);
    all = [all atril'];
 
    temp=[];
for t = 1:size(allcoh_2mWin,1);
        tt = 1:size(allcoh_2mWin,2);
    temp(t,tt,:)=squeeze(allcoh_2mWin(t,tt,1:821));
    all_sep = [all_sep squeeze(temp(t,tt,:))'];
end  
    
end
all_sep_no1s = all_sep(all_sep<1);
all_sep_no1s = reshape(all_sep_no1s,821,length(all_sep_no1s)/821);



figure; stdshade(all(10:end,:)',0.5,'k');set(gca, 'XScale', 'log')
figure; stdshade(all_sep(10:end,:)',0.5);set(gca,'XScale','log')
figure; 
for y = 1:size(all_sep,2);
    plot(all_sep(10:end,y));hold on;set(gca, 'XScale', 'log')
end

mean005to05 = mean(all_sep_no1s(10:83,:));
figure; histogram(mean005to05,80)

% figure; 
% histogram(mean(all
% test = mean(overallcoh_2mWin(1:20,:));
% for r = 1:length(test)
%     if test(r) < 0.4
%         overallcoh_2mWin(:,r) = NaN;
%     end
% end

% 
% cohplot = overallcoh_2mWin(~isnan(overallcoh_2mWin));
% cohplot=reshape(cohplot,100,length(cohplot)/100);
% figure; stdshade(cohplot',0.5,'k');set(gca, 'XScale', 'log')

%% Shuffled data control
% same axon
clear
prompt = {'Enter experiment type (ACh or NA or VCIN)'};
dlgtitle = 'Input';
dims = [1 35];

definput = {'ACh'};
exptype = inputdlg(prompt,dlgtitle,dims,definput);

Dir1coh = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Preprocessed2P\SmallFiles\SameAxon');
Dir2coh = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Analysis\WindowedCorrelation\SameAxon');


overallcoh = [];
overallcoh_shuf=[];


cd(Dir1coh)
SessCoh=uigetfile('*.mat','MultiSelect','on');
for x =1:length(SessCoh)
    cd(Dir1coh)
    clearvars -except SessCoh Dir1coh Dir2coh exptype x  overallcoh overallcoh_shuf
    load(char(SessCoh(x)),'Fnorm','filename');

r(1,:) = randperm(length(Fnorm));
r(2,:) = randperm(length(Fnorm));
Fnorm_shuf = Fnorm(r);

 Fnorm = lowpass(Fnorm',10,100);
 Fnorm_shuf = lowpass(Fnorm_shuf',10,100);
Fnorm=Fnorm';
Fnorm_shuf=Fnorm_shuf';

 f = linspace(0.01,5);
 win=18000;
 ov=win*0.98;
% [coh,fcoh] = mscohere(Fnorm(1,:),Fnorm(2,:),win,ov,f,100);
% [coh_shuf, fcoh]= mscohere(Fnorm_shuf(1,:),Fnorm_shuf(2,:),win,ov,f,100);

% below is to try axon-axon coherence like behXcorr
nfft = 12000; %2m window used to calculate fft
 noverlap =nfft*0.98; %98 percent overlap
 %[coh,fcoh] = mscohere(Fnorm(1,:),Fnorm(2,:),hann(nfft),noverlap);
%[coh_shuf, fcoh]= mscohere(Fnorm_shuf(1,:),Fnorm_shuf(2,:),hann(nfft),noverlap);
%below as of 5/4/22
 [coh_shuf,fcoh] = mscohere(Fnorm_shuf(1,:),Fnorm_shuf(2,:),hann(12000),12000*0.98,[],100);
 [coh,fcoh] = mscohere(Fnorm(1,:),Fnorm(2,:),hann(12000),12000*0.98,[],100);




cd(Dir2coh)

save(strcat('sameCoh_final20220504_LP10',filename(7:end-4)))

end
%%
cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Analysis\WindowedCorrelation\SameAxon'));
combCoh=uigetfile('sameCoh_*','MultiSelect','on');
allcoh=[];
allcoh_shuf=[];
for i = 1:length(combCoh)
    load(char(combCoh(i)))
allcoh=[allcoh coh];%transposing removed to run like BehXCoh
allcoh_shuf=[allcoh_shuf coh_shuf];
end

save('SameCoh_MultWin_likeBehXCorr_all')
% 
% figure; stdshade(allcoh(2:end,:)',0.5,'k');set(gca, 'XScale', 'log')
% figure; stdshade(allcoh_shuf(2:end,:)',0.5,'k');set(gca, 'XScale', 'log')

 figure; stdshade(allcoh(10:821,:)',0.5,'k');set(gca, 'XScale', 'log')
 hold on; stdshade(allcoh_shuf(10:821,:)',0.5,'k');set(gca, 'XScale', 'log')
 
 %% heterogeneity within coherence plot
 load('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\ACh\Analysis\WindowedCorrelation\difCoh_axonaxon_final20220504_10hzLP1693_190702_4roisplane3.mat')
allcoh=[];
%%
cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\ACh\Analysis\WindowedCorrelation\'));
Coh=uigetfile('difCoh_*','MultiSelect','on');
allcoh=[];

for i = 1:length(Coh)
    load(char(Coh(i)),'allcoh_2mWin')

allcoh_2mWin_tril=[]; 
for u = 1:821;
allcoh_2mWin_tril(:,:,u) = tril(squeeze(allcoh_2mWin(:,:,u)),-1);
end


for y = 1:size(allcoh_2mWin_tril,1)
    for yy = 1:size(allcoh_2mWin_tril,1)
        allcoh=[allcoh squeeze(allcoh_2mWin_tril(y,yy,:))];
    end
end
end

allcoh_unique = unique(allcoh','rows');

allcoh_use = allcoh_unique';

allcoh_m05to005 = mean(allcoh_use(10:82,:));

figure; histogram(allcoh_m05to005,80)

