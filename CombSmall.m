
function []= CombSmall(fn);
prompt = {'Enter experiment type (ACh or NA)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'ACh'};
Type = inputdlg(prompt,dlgtitle,dims,definput);

Fnorm_comb = [];
dFF_comb = [];
F_comb = [];
Fnorm_mean_comb=[];
NumAxons = [];
F_decon_comb=[];
F_decon_mean_comb=[];
cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(Type(1,1)),'\Preprocessed2P\SmallFiles'));
fn=uigetfile('*.mat','Select ALL FILES for a SNIGLE SESSION','MultiSelect','on'); % for multiple files

for i =  1:length(fn)% for multiple files
    %i=1; % for sinlge files
    load(char(fn(i))); % use this for multiple files
Fnorm_comb = [Fnorm_comb Fnorm'];
dFF_comb = [dFF_comb dFF'];
F_comb = [F_comb F'];
if size(Fnorm,1) ==1
    Fnorm_mean_comb = [Fnorm_mean_comb Fnorm];
    %F_decon_comb = [F_decon_comb F_decon];
end
if size(Fnorm,1) > 1
Fnorm_mean_comb = [Fnorm_mean_comb mean(Fnorm)];
%F_decon_mean_comb = [F_decon_mean_comb mean(F_decon)];
end
NumAxons=[NumAxons size(dFF,1)];
end

clearvars -except pupil whisk walk Fnorm_mean_comb F_comb dFF_comb Fnorm_comb filename X Y Z Ztop Type NumAxons F_decon_comb F_decon_mean_comb
num = round(length(Fnorm_mean_comb)/size(dFF_comb,1));
Fnorm_mean_comb = Fnorm_mean_comb(1:size(dFF_comb,1)*num);
Fnorm_mean_comb = reshape(Fnorm_mean_comb,size(dFF_comb,1),length(Fnorm_mean_comb)/size(dFF_comb,1));
%F_decon_mean_comb = F_decon_mean_comb(1:size(dFF_comb,1)*num);
%F_decon_mean_comb = reshape(F_decon_mean_comb,size(dFF_comb,1),length(F_decon_mean_comb)/size(dFF_comb,1));
cd(strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(Type(1,1)),'\Preprocessed2P\SmallFiles\Combined'));
save(filename(7:end-10));
end

