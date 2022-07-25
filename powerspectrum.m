%% power spectrum

clear

prompt = {'Enter experiment type (ACh or NA or VCIN)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'ACh'};
exptype = inputdlg(prompt,dlgtitle,dims,definput);

cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Preprocessed2P/SmallFiles/Combined'))
file = uigetfile('','MultiSelect','on');

ps_all = [];
    ps_whisk_all = [];
    ps_walk_all = [];
    ps_pupil_all = [];

    figure;
for r = 1:length(file)
    load(char(file(r)),'dFF_comb','whisk','walk','pupil')
  dFF_comb = dFF_comb(1:6000,:);
    dFF_comb=mean(dFF_comb,2);
    dFF_comb=highpass(dFF_comb,1,100);
    whisk = highpass(whisk(1:6000),1,100);
    walk = highpass(walk(1:6000),1,100);
    pupil = highpass(pupil(1:6000),1,100);
   
    [ps f] = pspectrum(dFF_comb);
    ps_whisk =pspectrum(whisk);
    ps_walk =pspectrum(walk);
    ps_pupil =pspectrum(pupil);
    %for t = 1:size(dFF_comb,2)
%         dFF_comb(1:6000,r) = highpass(dFF_comb(1:6000,r),1,100);
%         [ps(t,:) f] = pspectrum(dFF_comb(1:6000,t));
    %end
    
    ps_all = [ps_all ps];
    ps_whisk_all = [ps_whisk_all ps_whisk];
    ps_walk_all = [ps_walk_all ps_walk];
    ps_pupil_all = [ps_pupil_all ps_pupil];
    
    subplot(4,5,r); set(gca, 'XScale', 'log')
hold on;stdshade(ps_walk',0.5,'c')
stdshade(ps_whisk',0.5,'g')
stdshade(ps_pupil',0.5,'m')
end


f_corr = f*0.5*(100/3.14);

ps_all_corr = ps_all;
for u = 1:size(ps_all,2)
    ps_all_corr(:,u) = ps_all(:,u)/max(ps_all(:,u));
    ps_all_whisk_corr(:,u) = ps_whisk_all(:,u)/max(ps_whisk_all(:,u));
    ps_all_walk_corr(:,u) = ps_walk_all(:,u)/max(ps_walk_all(:,u));
    ps_all_pupil_corr(:,u) = ps_pupil_all(:,u)/max(ps_pupil_all(:,u));
end


figure; stdshade(ps_all_corr',0.5,'k')
hold on;
stdshade(ps_all_whisk_corr',0.5,'g')
stdshade(ps_all_walk_corr',0.5,'c')
stdshade(ps_all_pupil_corr',0.5,'m')
set(gca, 'XScale', 'log')
xlim([83,820])%0.5 to 10 hz

%% power spectrum for behavioral data

clear


cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\ACh\Preprocessed2P/SmallFiles/Combined'))
achfile = uigetfile('','MultiSelect','on');
cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\NA\Preprocessed2P/SmallFiles/Combined'))
nafile = uigetfile('','MultiSelect','on');
cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\VCIN\Preprocessed2P/SmallFiles/Combined'))
vcinfile = uigetfile('','MultiSelect','on');

file = horzcat(achfile,nafile,vcinfile);
 
    ps_whisk_all = [];
    ps_walk_all = [];
    ps_pupil_all = [];
    ps_whisk_control_all = [];
             ps_whisk_corr_all = [];
    ps_walk_corr_all = [];
    ps_pupil_corr_all = [];
    ps_whisk_control_corr_all = [];
%     psdx_whisk_all = [];
%     psdx_walk_all = [];
%     psdx_pupil_all = [];
figure;
for r = 1:length(file)
    load(char(file(r)),'whisk','walk','pupil')
    
    if length(whisk) > 6000
        whisk=highpass(whisk,1,100);
        walk=highpass(walk,1,100);
        pupil=highpass(pupil,1,100);
        
    t=0.01:0.01:60;
    S = 0.7*sin(2*pi*5*t) + sin(2*pi*.05*t);
    whisk_control=whisk(1:6000)+S';
    
        [ps_whisk f] = pspectrum(whisk(1:6000));
        [ps_walk f] = pspectrum(walk(1:6000));
        [ps_pupil f] = pspectrum(pupil(1:6000));
        [ps_whisk_control f] = pspectrum(whisk_control(1:6000));

    
    ps_whisk_all = [ps_whisk_all ps_whisk];
    ps_walk_all = [ps_walk_all ps_walk];
    ps_pupil_all = [ps_pupil_all ps_pupil];
    ps_whisk_control_all = [ps_whisk_control_all ps_whisk_control];
    
% x = whisk;
% N = length(x);
% xdft = fft(x);
% xdft = xdft(1:N/2+1);
% psdx_whisk = (1/(2*pi*N)) * abs(xdft).^2;
% psdx_whisk(2:end-1) = 2*psdx_whisk(2:end-1);
% freq_whisk = 0:(2*pi)/N:pi;
% 
% x = walk;
% N = length(x);
% xdft = fft(x);
% xdft = xdft(1:N/2+1);
% psdx_walk = (1/(2*pi*N)) * abs(xdft).^2;
% psdx_walk(2:end-1) = 2*psdx_walk(2:end-1);
% freq_walk = 0:(2*pi)/N:pi;
% 
% x = pupil;
% N = length(x);
% xdft = fft(x);
% xdft = xdft(1:N/2+1);
% psdx_pupil = (1/(2*pi*N)) * abs(xdft).^2;
% psdx_pupil(2:end-1) = 2*psdx_pupil(2:end-1);
% freq_pupil = 0:(2*pi)/N:pi;
% 
%     psdx_whisk_all = [psdx_whisk_all psdx_whisk(1:2000)];
%     psdx_walk_all = [psdx_walk_all psdx_walk(1:2000)];
%     psdx_pupil_all = [psdx_pupil_all psdx_pupil(1:2000)];

ps_whisk_corr = ps_whisk_all;

    ps_whisk_corr = ps_whisk_all/max(ps_whisk_all);
      ps_walk_corr = ps_walk_all/max(ps_walk_all);
        ps_pupil_corr = ps_pupil_all/max(ps_pupil_all);
        ps_whisk_control_corr = ps_whisk_control_all/max(ps_whisk_control_all);
        
          ps_whisk_corr_all = [ps_whisk_corr_all ps_whisk_corr];
    ps_walk_corr_all = [ps_walk_corr_all ps_walk_corr];
    ps_pupil_corr_all = [ps_pupil_corr_all ps_pupil_corr];
    ps_whisk_control_corr_all = [ps_whisk_control_corr_all ps_whisk_control_corr];


subplot(7,7,r); stdshade(ps_whisk_control_corr',0.5,'k');set(gca, 'XScale', 'log')
hold on;stdshade(ps_walk_corr',0.5,'c')
stdshade(ps_whisk_corr',0.5,'g')
stdshade(ps_pupil_corr',0.5,'m')
end
end
f_corr = f*0.5*(100/3.14);


figure; stdshade(ps_whisk_corr_all',0.5,'g');hold on
stdshade(ps_walk_corr_all',0.5,'c')
stdshade(ps_pupil_corr_all',0.5,'m')
stdshade(ps_whisk_control_corr_all',0.5,'k')
set(gca, 'XScale', 'log')
xlim([42,820])

% figure;
% plot(freq_whisk(1:2000)/pi,10*log10(mean(psdx_whisk_all,2)))
% hold on;plot(freq_walk(1:2000)/pi,10*log10(mean(psdx_walk_all,2)))
% plot(freq_pupil(1:2000)/pi,10*log10(mean(psdx_pupil_all,2)))
% grid on
% title('Periodogram Using FFT')
% xlabel('Normalized Frequency (\times\pi rad/sample)') 
% ylabel('Power/Frequency (dB/rad/sample)')
%% behavior this one's working?
clear


cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\ACh\Preprocessed2P/SmallFiles/Combined'))
achfile = uigetfile('','MultiSelect','on');
cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\NA\Preprocessed2P/SmallFiles/Combined'))
nafile = uigetfile('','MultiSelect','on');
cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\VCIN\Preprocessed2P/SmallFiles/Combined'))
vcinfile = uigetfile('','MultiSelect','on');

file = horzcat(achfile,nafile,vcinfile);
 
    ps_whisk_all = [];
    ps_walk_all = [];
    ps_pupil_all = [];
    psdx_whisk_all = [];
    psdx_walk_all = [];
    psdx_pupil_all = [];

for r = 1:length(file)
    load(char(file(r)),'whisk','walk','pupil')
    t=0.01:0.01:600;
    S = 0.7*sin(2*pi*5*t) + sin(2*pi*.05*t);
    whisk_60=whisk(1:60000)+S';
    if length(whisk) > 60000
    whisk_60 = whisk_60(1:60000);
    
Fs = 100;
X = whisk_60;
Y = fft(X);
L = length(whisk_60);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

P1_whisk(:,r) = timeseries(P1,f);

X = walk;
Y = fft(X);
L = length(whisk);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

P1_walk(:,r) = timeseries(P1,f);

X = pupil;
Y = fft(X);
L = length(whisk);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

P1_pupil(:,r) = timeseries(P1,f);

end

r = 0.01:0.01:100;
for i = 1:length(file)
 temp= resample(P1_whisk(i),r);
 P1_whisk_resam(:,i) = temp.data;
  temp= resample(P1_walk(i),r);
 P1_walk_resam(:,i) = temp.data;
  temp= resample(P1_pupil(i),r);
 P1_pupil_resam(:,i) = temp.data;
 f = temp.time;
end
end

figure; plot(f,mean(P1_whisk_resam,2),'g');
hold on;plot(f,mean(P1_walk_resam,2),'c');
plot(f,mean(P1_pupil_resam,2),'m');
set(gca, 'XScale', 'log')

%% now axons
clear

prompt = {'Enter experiment type (ACh or NA or VCIN)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'ACh'};
exptype = inputdlg(prompt,dlgtitle,dims,definput);

cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Preprocessed2P/SmallFiles/Combined'))
file = uigetfile('','MultiSelect','on');

P1_all = [];

for r = 1:length(file)
    load(char(file(r)),'dFF_comb')
 Fs=100;   
    for t = 1:size(dFF_comb,2)
        
       X = dFF_comb(:,t);
Y = fft(X);
L = length(dFF_comb);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
    
    

  P1 = timeseries(P1,f);
  tt =0.01:0.0001:10;
  temp= resample(P1,tt);
 P1_resam(:,t) = temp.data;
 f = temp.time;
 
 P1_all = [P1_all P1_resam];
end
end

figure; plot(f,median(P1_resam,2),'k');
set(gca, 'XScale', 'log')

%% try again with pspectrum
