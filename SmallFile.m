function []= SmallFile(Sessions);
prompt = {'Enter experiment type (ACh or NA or VCIN)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'ACh'};
exptype = inputdlg(prompt,dlgtitle,dims,definput);

Dir1small = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Preprocessed2P\Meso');
Dir2small = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Preprocessed2P\SmallFiles');
%  Dir1small = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Preprocessed2P\Meso\SameAxon');
%  Dir2small = strcat('\\ion-nas.uoregon.edu\mccormicklab2\Lindsay\ACh_NA_Synch_Project\',char(exptype(1,1)),'\Preprocessed2P\SmallFiles\SameAxon');

 
cd(Dir1small)
Sessions=uigetfile('*.mat','Select the INPUT DATA FILE(s)','MultiSelect','on');
for q =  1:length(Sessions)
    cd(Dir1small)
    load(char(Sessions(q)));
%Session=uigetfile('*.mat','Select Spike2 .mat file to preprocess');
%Session = cellstr(Session);

  %file = load(char(Session));
    %display('File loaded.')
    
walk = resam_LP_walk;
whisk =  resam_LP_whisk;
pupil = resam_LP_pupil;
Fnorm = resam_LP_Fnorm;
F = resam_LP_F;
dFF = resam_LP_dff;
X = meta_CenterX;
Y = meta_CenterY;
Z = meta_CenterZ;
Ztop = meta_Ztop;
%AdjustX = meta_AdjustX+X;
%AdjustY = meta_AdjustY+Y;



 this_Session = char(Sessions(q));
filename = strcat('Small_',this_Session(8:end));
 
clearvars -except AdjustX AdjustY Z X Y Ztop pupil whisk walk F dFF Fnorm filename Sessions Dir1small  Dir2small q
cd(Dir2small)
save(filename)


disp('Done.')
end
end