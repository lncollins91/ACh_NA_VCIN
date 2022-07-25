
%% ACh/NA project pipeline
% Analyze data and create figures for Collins, Francis, Emanuel, & McCormick 2022
% All code written by Lindsay Collins, updated 7/8/2022
%% Preprocessing
% Preprocess Spike2 file 
  % Can run as many files as you'd like at once
  % Saves a file in Preprocessed2P > Spike2
S2

% Preprocess 2P data
  % Be sure you've already selected ROIs in Suite2P!
  % Run each session separately 
  % Saves a file in Preprocessed2P > Meso
SNRDFF 

% Save SmallFile
  % Can run as many sessions as you'd like at once
SmallFile

% Combine SmallFiles
  % Run each session separately
CombSmall

% Whisk and Walk Onsets and Offsets
  % Can run all at once
WhiskOnOff
WalkOnOff

%% Cross correlations and coherence between axons/cells and behavioral data
  % Can run all at once, will get individual graphs
BehXCorr
BehXCorr_comb
BehXCoh

% Cross correlations between axons, sliding window
  % Can run all at once
WinCorr

% Cross correlations between axons as  function of arousal state
  % Can run all at once, writes over Corr file
%XCbyAr
XCbyBeh

% 
ShufCorr
ShufCorrWhisk 
CorrBoxplots %Fig. 4 B and E, need ShufCorr and ShufCorrWhisk to run
VCIN_pseudoraster % Fig. 6 F-H

%% Coherence between axons during walking and stillness
  % Can run multiple, takes a long time
wcohWS %not currently using

%% 
CorrByDist %this makes the last figures

%%
MovieGUI % this runs FaceMap
FaceMap1_working % this combines Facemap with dff, fnorm, pupil, walk, and whisk data

%% LINDSAY TO DO NEXT:
% Facemap data still isn't lining up correctly. wtf. 
% 
