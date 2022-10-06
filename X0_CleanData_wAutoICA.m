% master script
global DataConfig

% which experiment are we going to run?
ConfigFileName = 'Config_Danielle_051022';

% what do we need to do?
% IRRELEVANT NOW. With AutoICA it's all one mode.
ModeToPerform = 'AutoICA';
% 'PreICA' = preICA preparations, including decomposition and plot ICA
% 'PostICA' = remove ICA components, epoch and baseline

ResetICAcomponents = 1;
% default should be to overwrite.
% 1 = overwrite ICAcomponents.xlsx file when you're in preICA mode.
% 0 = do not overwrite ICAcomponents.xlsx.
% token change.

%% load analysis parameters
% first patch together full destination of the config file and add access
% to all the subfunctions you're about to call.
Current_File_Path = pwd;
addpath('Functions');
ConfigFilePath = [Current_File_Path filesep 'SupportingDocs' filesep ConfigFileName '.xlsx'];
Options = detectImportOptions(ConfigFilePath);

for k = 1:numel(Options.VariableTypes)
    Options.VariableTypes{k} = 'char';
end
DataConfig = table2struct(readtable(ConfigFilePath, Options));

% feed in info about what you want the config file to contain or do.
DataConfig.mode = ModeToPerform;
DataConfig.ResetICAcomponents = ResetICAcomponents;

% adjust the config data to suit our purposes. Essentially takes the data
% from the config file and makes a structure of cell arrays (mostly 1 x 1).
DataConfig = adjustConfigData(DataConfig);

% make standard AR parameter files (make this conditional on a "CustomAR"
% field of DataConfig in moment.

if DataConfig.CustomAR{1} == 0
    CreateIndividualARfiles(DataConfig);
else
    % need to create four files on your own.
    % ICA_Prep_Values.xlsx
    % AR_Parameters_for_MW_Blinks.xlsx
    % AR_Parameters_for_SVT_CRAP.xlsx
    % ICA_components.xlsx
end
%% The actual processing steps.
currentTime = clock;
% the 5th element of tVar is the current time in minutes.
% calculates minutes since midnight.
StartTime = currentTime(4)*60 + currentTime(5);
% 
% 
% %% do the preprocessing via PREP pipeline (or equivalent if PREP is
% % excluded). Also, can't do this parallel, as PREP calls up a
% % multithread loop, and those can't be nested.
% Y1_preprocess_wPREP;
% 
% %% prepare for next step
% DataConfig.LastProcess = cellstr('X1_PreProcess');
% DataConfig.LastSuffix = cellstr('_ds_addChans_PREP_bp_refs.set');
% 
% %% adjust events so they're on the correct timeline (if needed).
% tmpDataConfig = DataConfig;
% totalSUBS = length(tmpDataConfig.SUB);
% parfor loopIdx = 1:totalSUBS
%     SUB =  tmpDataConfig.SUB(loopIdx);
%     X1b_fixEvents_p(tmpDataConfig, SUB);
% end
% 
% %% prepare for next step
% DataConfig.LastProcess = cellstr('X1b_fixEvents');
% DataConfig.LastSuffix = cellstr('_ds_addChans_PREP_bp_refs_event.set');
% 
% %% do the ICA prep
% tmpDataConfig = DataConfig;
% totalSUBS = length(tmpDataConfig.SUB);
% parfor loopIdx = 1:totalSUBS
%     SUB =  tmpDataConfig.SUB(loopIdx);
%     X2_icaprep_p(tmpDataConfig, SUB);
% end
% 
% %% prepare for next step
% DataConfig.LastProcess = cellstr('X2_icaprep');
% DataConfig.LastSuffix = cellstr('_ds_addChans_PREP_bp_refs_event_icaPrep2.set');
% 
% %% and run the ICA decomp
% tmpDataConfig = DataConfig;
% totalSUBS = length(tmpDataConfig.SUB);
% parfor loopIdx = 1:totalSUBS
%     SUB =  tmpDataConfig.SUB(loopIdx);
%     X3_RunICA_p(tmpDataConfig, SUB);
% end
% 
% %% prepare for next step (can't update DataConfig in parallel).
% DataConfig.LastProcess = cellstr('X3_RunICA');
% DataConfig.LastSuffix = cellstr('_ds_addChans_PREP_bp_refs_event_icaWeighted.set');
% 
% %% plot the topos separately.
% % for some reason we can run out of memory if we try it together with
% % Runica.
% tmpDataConfig = DataConfig;
% totalSUBS = length(tmpDataConfig.SUB);
% for loopIdx = 1:totalSUBS
%     SUB =  tmpDataConfig.SUB(loopIdx);
%     X3b_PlotICAtopos_p(tmpDataConfig, SUB);
% end
% 
% %% prepare for next step (can't update DataConfig in parallel).
% DataConfig.LastProcess = cellstr('X3b_PlotICAtopos');
% DataConfig.LastSuffix = cellstr('_ds_addChans_PREP_bp_refs_event_icaWeighted.set');
% 
% %% remove the noisy components
% tmpDataConfig = DataConfig;
% totalSUBS = length(tmpDataConfig.SUB);
% for loopIdx = 1:totalSUBS
%     SUB =  tmpDataConfig.SUB(loopIdx);
%     X4_RemoveICA_p(tmpDataConfig, SUB);
% end
% 
% %% prepare for next step
% DataConfig.LastProcess = cellstr('X4_RemoveICA');
% DataConfig.LastSuffix = cellstr('_ds_PREP_ica_corr_cbip.set');
% 
% %% bin the epochs defined earlier.
% tmpDataConfig = DataConfig;
% totalSUBS = length(tmpDataConfig.SUB);
% parfor loopIdx = 1:totalSUBS
%     SUB =  tmpDataConfig.SUB(loopIdx);
%     X5_BinEpochs_p(tmpDataConfig, SUB);
% end
% 
% %% prepare for next step (can't update DataConfig in parallel).
% DataConfig.LastProcess = cellstr('X5_BinEpochs');
% DataConfig.LastSuffix = cellstr('_ds_PREP_ica_corr_cbip_elist_bins_epoch.set');
% 
% %% artifact rejection (according to config file).
% tmpDataConfig = DataConfig;
% totalSUBS = length(tmpDataConfig.SUB);
% for loopIdx = 1:totalSUBS
%     SUB =  tmpDataConfig.SUB(loopIdx);
%     imageType = 'none'; % or 'pdf' but this fails in parallel mode because it demands too much memory.
%     % and so does 'png' under enough load. So optionally 'none' too. 
%     %
%     X6_ArtifactRejection_p(tmpDataConfig, SUB, imageType);
% end

%% prepare for next step (can't update DataConfig in parallel).
DataConfig.LastProcess = cellstr('X6_ArtifactRejection');
DataConfig.LastSuffix = cellstr('_ds_PREP_ica_corr_cbip_elist_bins_epoch_ar.set');

%% and may as well extract the data here too.
tmpDataConfig = DataConfig;
totalSUBS = length(tmpDataConfig.SUB);
parfor loopIdx = 1:totalSUBS
    SUB =  tmpDataConfig.SUB(loopIdx);
    X7_ExtractEpochedData_p(tmpDataConfig, SUB);
end

% clock yields 6 element array of date/time info.
% element 4 is time in hours, 5 is minutes.
currentTime = clock;
% calculates minutes since midday/midnight.
EndTime = currentTime(4)*60 + currentTime(5);
disp(['Time taken for total analysis ' num2str((EndTime-StartTime)/60) ' minutes']);