% ANALYSIS PIPELINE FOR EEG-DATA - PREPROCESSING
% Epoch data (continuous data into trials)
% Last edited 20-06-2016

function eeg_run_trial_definition_submean(sub_nr)
% Input:    - sub_nr:           subject number
%           - EEG-data (preprocessed)
% 
% Output:   - sb_??_eegdata_preprocessed_trials.mat
%           - trial_matrix.mat


%% DIRECTORIES & FILES
% General
addpath('/home/valmon/SUBMEAN/Scripts_analyses/my_scripts_group3');
dirs.home           = 'U:\SUBMEAN\';

% dirs.home           = '/home/valmon/SUBMEAN/';

% Subject specific folders
dirs.sb.data        = fullfile(dirs.home,'SUBMEAN_Data_group3','EEG','data',['sb_',num2str(sub_nr)]);
dirs.sb.save        = dirs.sb.data;
dirs.sb.info        = fullfile(dirs.sb.save,'info');

dirs.sb.triggers =     fullfile(dirs.home,'SUBMEAN_Data_group3','EEG','Raw_data',['sb_',num2str(sub_nr)]);

% Subject filenames
files.sb.data       = ['sb_',num2str(sub_nr),'_eegdata_preprocessed.mat'];
files.sb.save       = ['sb_',num2str(sub_nr),'_eegdata_preprocessed_trials.mat'];
files.sb.triggers   = dir(fullfile(dirs.sb.triggers,['*.vmrk']));
files.sb.triggers   = [files.sb.triggers.name];


%% LOAD DATA
% load preprocessed EEG data
load(fullfile(dirs.sb.data,files.sb.data)); 
% Load eeg marker and their time stamp from trigger logfile (.vmrk)
triggers=bva_readmarker(fullfile(dirs.sb.triggers,files.sb.triggers)); 
% Transpose the matrix
triggers=triggers';

% Define relevant triggers
picture_triggers=[30:45,106:109,130:145]; 

% Store the row numbers corresponding to the relevant triggers
picture_list_ids=[];
for t=1:length(picture_triggers)
    picture_list_ids=[picture_list_ids;find(triggers(:,1)==picture_triggers(t))];
end


%% COMPILE TRL MATRIX 

% Compile trial matrix 1
% 1. Time stamp of the trigger
% 2. Block type portcode
% 3. Visibility portcode
% 4. Blank 1 or mask 1 portcode
% 5. Prime 1 portcode
% 6. Blank 2 or mask 2 portcode
% 7. Prime 2 portcode
% 8. Blank 3 or mask 3 portcode
% 9. Prime 3 portcode
% 10. Blank 4 or mask 4 portcode (or position 1/2/3 portcode for noun or verb codes)
% 11. Picture portcode
% 12. Congruency portcode
% 13. Response portcode
% 14. Trial number

trial_matrix=[];
for t=1:length(picture_list_ids)
    trial_matrix(t,1)=triggers(picture_list_ids(t),2);
    trial_matrix(t,2:13)=triggers(picture_list_ids(t)-9:picture_list_ids(t)+2,1);
end

% Sort the matrix by the time stamp of the trigger
[~,i]=sort(trial_matrix(:,1));
trial_matrix=trial_matrix(i,:);

% Add trial number column
trial_matrix(:,14)=1:length(trial_matrix);

% Check whether the matrix has the right number of trials
if  size(trial_matrix,1) ~= 1056
    error('Problem with amount of trials!');
end

% Delete trials with timing problem
trials_timingproblem = delete_timingproblems(sub_nr);
trial_matrix(trials_timingproblem,:)=[];

% Throw away trials with response problem
trial_matrix((trial_matrix(:,13)<3) | (trial_matrix(:,13)>4),:)=[];

% Compile trial matrix 2
% 1. Onset sample
% 2. End sample
% 3. Trigger offset (negative value = onset before trigger)
% 4. Block type code (1: sentence; 2: noun; 3: verb)
% 5. Visibility code(1: visible; 0: invisible)
% 6. Congruency code (1: INcongruent; 0: congruent)
% 7. Response (1=match, 2=no match)
% 8. Correct/incorrect answer column (correct=1)
% 9. Trial number

% Time-lock to picture onset
pretrig_dur     = 52/60;     % pre-trigger is trial onset to picture onset           
posttrig_dur    = 1;         % post trigger is 1s = picture duration

% Create trl matrix
trl = nan(length(trial_matrix),1);

% Compute timing (columns 1-3)
trl(:,1) = trial_matrix(:,1)-ceil(pretrig_dur*data.fsample);
trl(:,2) = trial_matrix(:,1)+ceil(posttrig_dur*data.fsample);
trl(:,3) = -ceil(pretrig_dur*data.fsample);

% Add simpler condition codes
% Block type code (column 4)
trl(trial_matrix(:,2)==97,4)=1;
trl(trial_matrix(:,2)==98,4)=2;
trl(trial_matrix(:,2)==99,4)=3;
% Visibility codes (visible=1) (column 5)
% NB: visibility portcode error for verb blocks is fixed here
trl(trial_matrix(:,3)==55 | trial_matrix(:,3)==115 | trial_matrix(:,3)==151,5)=0;
trl(trial_matrix(:,3)==54 | trial_matrix(:,3)==114 | ((trial_matrix(:,2)==99)&(trial_matrix(:,4)==201)),5)=1;
% Congruency codes (INcongruent=1) (column 6)
trl(trial_matrix(:,12)==51 | trial_matrix(:,12)==52 | trial_matrix(:,12)==53 | trial_matrix(:,12)==117 | trial_matrix(:,12)==153,6)=1;
trl(trial_matrix(:,12)==50 | trial_matrix(:,12)==116 | trial_matrix(:,12)==152,6)=0;

% Response (1=match, 2=no match)(7)
trl(:,7)=trial_matrix(:,13)-2;
trl(trl(:,7)<=0,7)=NaN;

% Add correct/incorrect answer column (correct=1)
trl((trl(:,7)==1 & trl(:,6)==0) | (trl(:,7)==2 & trl(:,6)==1),8)=1;
trl((trl(:,7)==2 & trl(:,6)==0) | (trl(:,7)==1 & trl(:,6)==1),8)=0;
trl(isnan(trl(:,7)),8)=NaN;

% Add trial number
trl(:,9)=trial_matrix(:,14);

% Add original matrix
trl(:,10:21)=trial_matrix(:,2:13);

% Save matrix
save(fullfile(dirs.sb.info,'trial_matrix.mat'),'trl');


%% REDEFINE TRIALS
% Redefine
cfg             = [];
cfg.trl         = trl;
data            = ft_redefinetrial(cfg,data);



%% SAVE EPOCHED DATA
save(fullfile(dirs.sb.save,files.sb.save),'data'); %tetest;