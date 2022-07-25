
%% Whole session

clear

prompt = {'Enter experiment type (ACh or NA or VCIN)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'ACh'};
exptype = inputdlg(prompt,dlgtitle,dims,definput);

% Different axon data
cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Preprocessed2P\SmallFiles\Combined'))
pcfile = uigetfile('*.mat','MultiSelect','on');

for y = 1:length(pcfile)
    
    cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Preprocessed2P\SmallFiles\Combined'))
    load(char(pcfile(y)),'dFF_comb','filename','whisk','walk')
    
    cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Analysis\WindowedCorrelation'))
    load(strcat('CxD_axonaxon_',filename(7:end-10)),'corrs')

% figure;plot(whisk);ylim([0,0.1])
% [xcut ycut] = ginput();
% close all

whisk_corr = (whisk-min(whisk))/(max(whisk)-min(whisk));
whisk_low = whisk_corr;
whisk_cut = median(whisk_corr);
whisk_low(whisk_low>whisk_cut) = NaN;

dFF_nowhisk=dFF_comb;
whisk_nowhisk=whisk;
whisk_nowalk=whisk;
whisk_whole=whisk;

for n = 1:size(dFF_comb,1)
    if whisk_corr(n) >whisk_cut;
        dFF_nowhisk(n,:) = NaN;
        whisk_nowhisk(n)=NaN;
    end
end

dFF_nowhisk=dFF_nowhisk(~isnan(dFF_nowhisk));
dFF_nowhisk=reshape(dFF_nowhisk,length(dFF_nowhisk)/size(dFF_comb,2),size(dFF_comb,2));
whisk_nowhisk=whisk_nowhisk(~isnan(whisk_nowhisk));

  walk_corr=walk-median(walk);
walk_low = walk_corr;
walk_cut = 0.025;
walk_low(walk_low>walk_cut) = NaN;

dFF_nowalk=dFF_comb;

for n = 1:size(dFF_comb,1)
    if walk_corr(n) >walk_cut;
        dFF_nowalk(n,:) = NaN;
        whisk_nowalk(n) = NaN;
    end
end

dFF_nowalk=dFF_nowalk(~isnan(dFF_nowalk));
dFF_nowalk=reshape(dFF_nowalk,length(dFF_nowalk)/size(dFF_comb,2),size(dFF_comb,2));
whisk_nowalk=whisk_nowalk(~isnan(whisk_nowalk));
   
  
    
    for a = 1:size(dFF_nowhisk,2)
    for b = 1:size(dFF_nowhisk,2)
        dFF_nowhisk_temp= dFF_nowhisk;
        dFF_nowhisk_temp(:,a) = NaN;
        dFF_nowhisk_temp(:,b) = NaN;
        z = nanmean(dFF_nowhisk_temp,2);
        z = whisk_nowhisk;
        [PartCorr_func_nowhisk(a,b) pval_nowhisk(a,b)] = partialcorr(dFF_nowhisk(:,a),dFF_nowhisk(:,b),z);
        temp = corrcoef(dFF_nowhisk(:,a),dFF_nowhisk(:,b));
        corrs_nowhisk(a,b) = temp(1,2);
    end
    end
    

        for a = 1:size(dFF_nowalk,2)
    for b = 1:size(dFF_nowalk,2)
        dFF_nowalk_temp= dFF_nowalk;
        dFF_nowalk_temp(:,a) = NaN;
        dFF_nowalk_temp(:,b) = NaN;
        %z = nanmean(dFF_nowalk_temp,2);
        z=whisk_nowalk;
        [PartCorr_func_nowalk(a,b) pval_nowalk(a,b)] = partialcorr(dFF_nowalk(:,a),dFF_nowalk(:,b),z);
        temp = corrcoef(dFF_nowalk(:,a),dFF_nowalk(:,b));
        corrs_nowalk(a,b) = temp(1,2);
    end
        end
  
        for a = 1:size(dFF_comb,2)
    for b = 1:size(dFF_nowhisk,2)
        dFF_whole_temp= dFF_comb;
        dFF_whole_temp(:,a) = NaN;
        dFF_whole_temp(:,b) = NaN;
        %z = nanmean(dFF_whole_temp,2);
        z=whisk_whole;
        [PartCorr_func_whole(a,b) pval_whole(a,b)] = partialcorr(dFF_comb(:,a),dFF_comb(:,b),z);
        temp = corrcoef(dFF_comb(:,a),dFF_comb(:,b));
        corrs_whole(a,b) = temp(1,2);
    end
    end
    
    cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Analysis\PartialCorrelation'))
    save(strcat('PartCorrs_comb_whiskcommon',filename(7:end-10)))
    clearvars -except exptype y pcfile
end

%% same axon
% same axon with dFF and whisking as common signal

clear

prompt = {'Enter experiment type (ACh or NA or VCIN)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'ACh'};
exptype = inputdlg(prompt,dlgtitle,dims,definput);

% same axon data
cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Preprocessed2P\SmallFiles\SameAxon'))
pcfile = uigetfile('*.mat','MultiSelect','on');

 PartCorr_func_same_all = [];
        pval_same_all=[];
        corr_same_all=[];
        mean_whisk_all=[];
        mean_walk_all=[];
        mean_pupil_all=[];
for y = 1:length(pcfile)
    
    clearvars -except exptype y PartCorr_func_same_all pval_same_all corr_same_all pcfile mean_whisk_all mean_walk_all mean_pupil_all
    cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Preprocessed2P\SmallFiles\SameAxon'))
    load(char(pcfile(y)),'dFF','filename','whisk','walk','pupil')
    
    cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Preprocessed2P\SmallFiles\Combined'))
    load(strcat(filename(7:end-10),'.mat'),'dFF_comb')
    dFF_all = dFF_comb; % remove unless doing vcin
    %z = mean(dFF_all,2);
    z=whisk;
     
    
        [PartCorr_func pval] = partialcorr2(dFF(1,:)',dFF(2,:)',z);
        corr = corrcoef(dFF(1,:),dFF(2,:));
        corr=corr(1,2);
        
        mean_whisk = mean(whisk);
        mean_walk = mean(walk);
        mean_pupil = mean(pupil);
        
        PartCorr_func_same_all = [PartCorr_func_same_all PartCorr_func];
        pval_same_all=[pval_same_all pval];
        corr_same_all=[corr_same_all corr];
        mean_whisk_all = [mean_whisk_all mean_whisk];
        mean_walk_all=[mean_walk_all mean_walk];
        mean_pupil_all = [mean_pupil_all mean_pupil];
  

end

figure; scatter(corr_same_all,PartCorr_func_same_all,'o','k');xlim([-1,1]);ylim([-1,1]);hold on;line([-1,1],[-1,1]);xlabel('Ordinary Correlation');ylabel('Partial Correlation');

cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Analysis\PartialCorrelation'))
load('PartCorrs_all')
p_bi = All_p;
p_bi(All_p <0.01) = 1;
p_bi(All_p > 0.01) = 0;

plot_part = All_part(p_bi ==1);
plot_ord = All_ord(p_bi==1);

figure; boxplot(All_part);ylim([-0.6 1]);figure; boxplot(PartCorr_func_same_all);ylim([-0.6 1]);
[t,p,h,c]=ttest2(plot_part,PartCorr_func_same_all);

figure; 
subplot(3,1,1);scatter(PartCorr_func_same_all,mean_whisk_all);
subplot(3,1,2);scatter(PartCorr_func_same_all,mean_walk_all);
subplot(3,1,3);scatter(PartCorr_func_same_all,mean_pupil_all);
save('PartCorrs_all_sameaxon')
%% combine and plot for figure
clear
%select files in this order: nowalk, whole, nowhisk
prompt = {'Enter experiment type (ACh or NA or VCIN)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'ACh'};
exptype = inputdlg(prompt,dlgtitle,dims,definput);

cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Analysis\PartialCorrelation'))
combfile = uigetfile('PartCorrs_comb','MultiSelect','on');
    All_part_whole =[];
    All_ord_whole = [];
    All_p_whole=[]
    All_part_mean_whole=[];
    All_ord_mean_whole=[];
    
        All_part_nowhisk =[];
    All_ord_nowhisk = [];
    All_p_nowhisk=[]
    All_part_mean_nowhisk=[];
    All_ord_mean_nowhisk=[];
    
        All_part_nowalk =[];
    All_ord_nowalk = [];
    All_p_nowalk=[]
    All_part_mean_nowalk=[];
    All_ord_mean_nowalk=[];
    
for g = 1:length(combfile)
%     load(char(combfile(g)),'PartCorr_func','corrs','pval','whisk')
%     PartCorr_func = reshape(PartCorr_func,1,size(PartCorr_func,1)*size(PartCorr_func,2));
%     corrs = reshape(corrs,1,size(corrs,1)*size(corrs,2));
%     pval=reshape(pval,1,size(pval,1)*size(pval,2));
%     All_part = [All_part PartCorr_func];
%     All_part_mean=[All_part_mean mean(PartCorr_func)];
%     All_ord_mean=[All_ord_mean mean(corrs)];
%     All_ord = [All_ord corrs];
%     All_p=[All_p pval];
%     %All_walk=[All_walk repmat(mean(walk),length(pval),1)'];
%     
%     clearvars PartCorr_func corrs pval whisk
% end
% cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Analysis\PartialCorrelation'))
%     save('PartCorrs_All_nowalk')
    
%% now pull in session mean for whole session
% cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Analysis\PartialCorrelation'))
% disp('Whole,nowhisk, nowalk')
% combfile = uigetfile('*.mat','MultiSelect','on');
%     All_part_whole =[];
%     All_ord_whole = [];
%     All_p_whole=[];
%     All_walk_whole=[];
%     All_part_mean_whole=[];
%     All_ord_mean_whole=[];
% for g = 1:length(combfile)
    load(char(combfile(g)),'PartCorr_func_whole','PartCorr_func_nowalk','PartCorr_func_nowhisk','corrs_whole','corrs_nowhisk','corrs_nowalk','pval_whole','pval_nowhisk','pval_nowalk')
   
    PartCorr_func_whole = reshape(PartCorr_func_whole,1,size(PartCorr_func_whole,1)*size(PartCorr_func_whole,2));
    corrs_whole = reshape(corrs_whole,1,size(corrs_whole,1)*size(corrs_whole,2));
    pval_whole=reshape(pval_whole,1,size(pval_whole,1)*size(pval_whole,2));
    All_part_whole = [All_part_whole PartCorr_func_whole];
    All_part_mean_whole=[All_part_mean_whole nanmean(PartCorr_func_whole)];
    All_ord_mean_whole=[All_ord_mean_whole mean(corrs_whole)];
    All_ord_whole = [All_ord_whole corrs_whole];
    All_p_whole=[All_p_whole pval_whole];
    %All_walk_whole=[All_walk_whole repmat(mean(walk),length(pval),1)'];
      PartCorr_func_nowhisk = reshape(PartCorr_func_nowhisk,1,size(PartCorr_func_nowhisk,1)*size(PartCorr_func_nowhisk,2));
    corrs_nowhisk = reshape(corrs_nowhisk,1,size(corrs_nowhisk,1)*size(corrs_nowhisk,2));
    pval_nowhisk=reshape(pval_nowhisk,1,size(pval_nowhisk,1)*size(pval_nowhisk,2));
    All_part_nowhisk = [All_part_nowhisk PartCorr_func_nowhisk];
    All_part_mean_nowhisk=[All_part_mean_nowhisk nanmean(PartCorr_func_nowhisk)];
    All_ord_mean_nowhisk=[All_ord_mean_nowhisk mean(corrs_nowhisk)];
    All_ord_nowhisk = [All_ord_nowhisk corrs_nowhisk];
    All_p_nowhisk=[All_p_nowhisk pval_nowhisk];
    
        PartCorr_func_nowalk = reshape(PartCorr_func_nowalk,1,size(PartCorr_func_nowalk,1)*size(PartCorr_func_nowalk,2));
    corrs_nowalk = reshape(corrs_nowalk,1,size(corrs_nowalk,1)*size(corrs_nowalk,2));
    pval_nowalk=reshape(pval_nowalk,1,size(pval_nowalk,1)*size(pval_nowalk,2));
    All_part_nowalk = [All_part_nowalk PartCorr_func_nowalk];
    All_part_mean_nowalk=[All_part_mean_nowalk nanmean(PartCorr_func_nowalk)];
    All_ord_mean_nowalk=[All_ord_mean_nowalk mean(corrs_nowalk)];
    All_ord_nowalk = [All_ord_nowalk corrs_nowalk];
    All_p_nowalk=[All_p_nowalk pval_nowalk];
    
    clearvars PartCorr_func_whole corrs_whole pval_whole PartCorr_func_nowhisk corrs_nowhisk pval_nowhisk PartCorr_func_nowalk corrs_nowalk pval_nowalk
end

% %now pull in esssion mean for nowhisk
% cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Analysis\PartialCorrelation'))
% combfile = uigetfile('PartCorrs_func_nowhisk*','MultiSelect','on');
%     All_part_nowhisk =[];
%     All_ord_nowhisk = [];
%     All_p_nowhisk=[];
%     All_walk_nowhisk=[];
%     All_part_mean_nowhisk=[];
%     All_ord_mean_nowhisk=[];
% for g = 1:length(combfile)
%     load(char(combfile(g)),'PartCorr_func','corrs','pval','whisk')
%     PartCorr_func = reshape(PartCorr_func,1,size(PartCorr_func,1)*size(PartCorr_func,2));
%     corrs = reshape(corrs,1,size(corrs,1)*size(corrs,2));
%     pval=reshape(pval,1,size(pval,1)*size(pval,2));
%     All_part_nowhisk = [All_part_nowhisk PartCorr_func];
%     All_part_mean_nowhisk=[All_part_mean_nowhisk mean(PartCorr_func)];
%     All_ord_mean_nowhisk=[All_ord_mean_nowhisk mean(corrs)];
%     All_ord_nowhisk = [All_ord_nowhisk corrs];
%     All_p_nowhisk=[All_p_nowhisk pval];
%     %All_walk_nowhisk=[All_walk_nowhisk repmat(mean(walk),length(pval),1)'];
%     
%     clearvars PartCorr_func corrs pval whisk
% end
% 
% 
% cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Analysis\PartialCorrelation'))
% combfile = uigetfile('PartCorrs_func_nowalk*','MultiSelect','on');
%     All_part_nowalk =[];
%     All_ord_nowalk = [];
%     All_p_nowalk=[];
%     %All_walk_nowalk=[];
%     All_part_mean_nowalk=[];
%     All_ord_mean_nowalk=[];
% for g = 1:length(combfile)
%     load(char(combfile(g)),'PartCorr_func','corrs','pval','whisk')
%     PartCorr_func = reshape(PartCorr_func,1,size(PartCorr_func,1)*size(PartCorr_func,2));
%     corrs = reshape(corrs,1,size(corrs,1)*size(corrs,2));
%     pval=reshape(pval,1,size(pval,1)*size(pval,2));
%     All_part_nowalk = [All_part_nowalk PartCorr_func];
%     All_part_mean_nowalk=[All_part_mean_nowalk mean(PartCorr_func)];
%     All_ord_mean_nowalk=[All_ord_mean_nowalk mean(corrs)];
%     All_ord_nowalk = [All_ord_nowalk corrs];
%     All_p_nowalk=[All_p_nowalk pval];
%     %All_walk_nowalk=[All_walk_nowalk repmat(mean(walk),length(pval),1)'];
%     
%     clearvars PartCorr_func corrs pval whisk
% end
% 


clearvars -except All*

diff_partOrd_nowalk = All_part_mean_nowalk-All_ord_mean_nowalk;
diff_partOrd = All_part_mean_whole-All_ord_mean_whole;
diff_partOrd_nowhisk = All_part_mean_nowhisk-All_ord_mean_nowhisk;

data = [diff_partOrd' diff_partOrd_nowalk' diff_partOrd_nowhisk']; %// Create row vector with your data
group = {'whole','nowalk','nowhisk'}; %// set the groups according to the data above

[p,tbl,stats] = anova1(data, group);
multcompare(stats)


figure; 
subplot(3,1,1);bar(mean(diff_partOrd));hold on;errorbar(mean(diff_partOrd),std(diff_partOrd));title('nowalk')
hold on;scatter(ones(length(diff_partOrd),1),diff_partOrd);ylim([-1,0])
subplot(3,1,2);bar(mean(diff_partOrd_nowalk));hold on;errorbar(mean(diff_partOrd_nowalk),std(diff_partOrd_nowalk));title('whole')
hold on;scatter(ones(length(diff_partOrd_nowalk),1),diff_partOrd_nowalk);ylim([-1,0])
subplot(3,1,3);bar(mean(diff_partOrd_nowhisk));hold on;errorbar(mean(diff_partOrd_nowhisk),std(diff_partOrd_nowhisk));title('nowhisk')
hold on;scatter(ones(length(diff_partOrd_nowhisk),1),diff_partOrd_nowhisk);ylim([-1,0])

figure; scatter(All_ord_whole,All_part_whole)

%% partial correlation x distance
clear
prompt = {'Enter experiment type (ACh or NA or VCIN)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'ACh'};
exptype = inputdlg(prompt,dlgtitle,dims,definput);
cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype),'\Analysis\WindowedCorrelation'));
pxdfile = uigetfile('CxD_axonaxon*','MultiSelect','on');

PartCorr_all = [];
corrs_all=[];
dists_all=[];

for i = 1:length(pxdfile)
    load(char(pxdfile(i)),'dFF_comb','corrs_tril','dists_tril');
    
    for a = 1:size(dFF_comb,2)
    for b = 1:size(dFF_comb,2)
        dFF_comb_temp= dFF_comb;
        dFF_comb_temp(:,a) = NaN;
        dFF_comb_temp(:,b) = NaN;
        z = nanmean(dFF_comb_temp,2);
        [PartCorr(a,b) pval(a,b)] = partialcorr(dFF_comb(:,a),dFF_comb(:,b),z);
    end
    end
    
    PartCorr_tril = tril(PartCorr,-1);
    PartCorr_tril(PartCorr_tril==0) = NaN;
    PartCorr_tril=PartCorr_tril(~isnan(PartCorr_tril));
    
    if length(PartCorr_tril) == length(corrs_tril)
    PartCorr_all = [PartCorr_all PartCorr_tril'];
    corrs_all = [corrs_all corrs_tril'];
    dists_all = [dists_all dists_tril'];
    end
    
end

mean_PartCorr_sameROI = nanmean(PartCorr_all(dists_all < 1));
mean_PartCorr_diffROI = nanmean(PartCorr_all(dists_all > 100));

std_PartCorr_sameROI = std(PartCorr_all(dists_all < 1));
std_PartCorr_diffROI = std(PartCorr_all(dists_all > 100));


save('PartCorrXDist_all')
    
 figure; scatter(dists_all,PartCorr_all,'.')  
 
 sameroi_mean = mean(PartCorr_all(dists_all < 1));
  sameroi_std = std(PartCorr_all(dists_all < 1));
 closeroi_mean = mean(PartCorr_all(dists_all < 1000 & dists_all > 1));
  closeroi_std = std(PartCorr_all(dists_all < 1000 & dists_all > 1));
 farroi_mean = mean(PartCorr_all(dists_all > 1000));
  farroi_std = std(PartCorr_all(dists_all > 1000));
   
    
%% partial correlation x size of signal (dF/F range)
clear
prompt = {'Enter experiment type (ACh or NA or VCIN)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'ACh'};
exptype = inputdlg(prompt,dlgtitle,dims,definput);
cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype),'\Analysis\WindowedCorrelation'));
pxdfile = uigetfile('CxD_axonaxon*','MultiSelect','on');

PartCorr_all = [];
OrdCorr_all=[];
Signal_range_all=[];
Signal_std_all=[];

for i = 1:length(pxdfile)
    load(char(pxdfile(i)),'dFF_comb');
    
    for a = 1:size(dFF_comb,2)
    for b = 1:size(dFF_comb,2)
        dFF_comb_temp= dFF_comb;
        dFF_comb_temp(:,a) = NaN;
        dFF_comb_temp(:,b) = NaN;
        z = nanmean(dFF_comb_temp,2);
        [PartCorr(a,b) pval(a,b)] = partialcorr(dFF_comb(:,a),dFF_comb(:,b),z);
        temp = corrcoef(dFF_comb(:,a),dFF_comb(:,b));
        OrdCorr(a,b) = corrcoef(1,2);
        rangea = max(dFF_comb(:,a))-median(dFF_comb(:,a));
        rangeb = max(dFF_comb(:,b))-median(dFF_comb(:,b));
        Signal_range(a,b) = mean([rangea,rangeb]);
        rangea = std(dFF_comb(:,a));
        rangeb = std(dFF_comb(:,b));
        Signal_std(a,b) = mean([rangea,rangeb]);
    end
    end
    
    PartCorr_tril = tril(PartCorr,-1);
    PartCorr_tril(PartCorr_tril==0) = NaN;
    PartCorr_tril=PartCorr_tril(~isnan(PartCorr_tril));
    
    Signal_range_tril = tril(Signal_range,-1);
    Signal_range_tril(Signal_range_tril==0) = NaN;
    Signal_range_tril=Signal_range_tril(~isnan(Signal_range_tril));
    
    Signal_std_tril = tril(Signal_std,-1);
    Signal_std_tril(Signal_std_tril==0) = NaN;
    Signal_std_tril=Signal_std_tril(~isnan(Signal_std_tril));
    
    OrdCorr_tril = tril(OrdCorr,-1);
    OrdCorr_tril(OrdCorr_tril==0) = NaN;
    OrdCorr_tril=OrdCorr_tril(~isnan(OrdCorr_tril));
    
    
    
   
    PartCorr_all = [PartCorr_all PartCorr_tril'];
    OrdCorr_all = [OrdCorr_all OrdCorr_tril'];
    Signal_std_all = [Signal_std_all Signal_std_tril'];
    Signal_range_all = [Signal_range_all Signal_range_tril'];
    
end

save('PartCorrXSignal_all')
    
 figure; scatter(Signal_range_all,PartCorr_all,'.');title('Range')
  figure; scatter(Signal_std_all,PartCorr_all,'.');title('std')