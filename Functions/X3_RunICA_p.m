%Script #3
%Operates on individual subject data
%Uses the output from Script #2: ICA_Prep.m
%This script loads the outputted semi-continuous EEG data file from Script #2, computes the ICA weights that will be used for artifact correction of ocular artifacts, transfers the ICA weights to the
%continuous EEG data file outputted from Script #1 (e.g., without the break periods and noisy segments of EEG removed), and saves a pdf of the topographic maps of the ICA weights.


% PLEASE NOTE:
% The results of ICA decomposition using binica/runica (i.e., the ordering of the components, the scalp topographies, and the time courses of the components) will differ slightly each time ICA weights are computed.
% This is because ICA decomposition starts with a random weight matrix (and randomly shuffles the data order in each training step), so the convergence is slightly different every time it is run.
% As a result, the topographic maps of the ICA weights and the excel spreadsheet (ICA_Components_MMN.xlsx) containing the list of ICA components to be removed for each subject included in this package
% will NOT be valid if ICA weights are re-computed. To avoid confusion or accidental overwriting of the relevant data files, this script has been commented out.

% To maintain the component weights and ordering from the original analysis, you can skip running this script and proceed to Script #4 Remove_ICA_Components.m.

% If you wish to re-compute ICA weights on the ERP CORE data, you will need to disregard the information in ICA_Components_MMN.xslx and evaluate the scalp topography and
% time course of the outputted ICA components to determine which component(s) to remove.

% To use this script, select all and use the shortcut Ctrl-T for PC or Command-T for Mac to uncomment the code.

function X3_AutoRunTheICA_p(DataConfig,SUB)

% initialize
close all;
clearvars -except DataConfig SUB;

    % Location of the main study directory
    DIR = fileparts(pwd);
    
    % location of preprocessing files.
    Current_File_Path = pwd;
    
    %***********************************************************************************************************************************************
    
    %Loop through each subject listed in SUB
    for i = 1:length(SUB)
        
        %Open EEGLAB and ERPLAB Toolboxes
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        
        %Define subject path based on study directory and subject ID of current subject
        Subject_Path = [DIR filesep SUB{i} filesep];
        
        % record where we're up to in case of crash.
        DataConfig.CurrentSUB = SUB{i};
        
        % Load the semi-continuous EEG data file outputted from Script #2 in .set EEGLAB file format
        FileToOpen = [SUB{i} DataConfig.LastSuffix{1}];
        EEG = pop_loadset( 'filename', FileToOpen, 'filepath', Subject_Path);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [SUB{i} '_ds_addChans_PREP_bp_refs_event_icaPrep2'], 'gui', 'off');
        
        % Compute ICA weights with runICA.
        % Try to include reference channels, so if mastoid referenced, then
        % add them to the ICA. (if not, exclude all non-scalp). 
        if strcmp(DataConfig.ReReference{1}, 'Mastoid') 
                    ChansForICA = [DataConfig.firstScalp:DataConfig.lastScalp, ...
            DataConfig.KeyChans{3}, DataConfig.KeyChans{4}]; 
        else 
            % no need to include mastoids if they aren't used for referencing
            ChansForICA = [DataConfig.firstScalp:DataConfig.lastScalp];
        end

        % runica can sometimes miscalculate rank and screw up. So make it
        % guess the rank, and force that rank onto the calculation using
        % 'pca' key value pair to avoid complex components that break things.
        EffectiveRank = rank(EEG.data(ChansForICA,:));
        
        EEG = pop_runica(EEG,'extended',1,'chanind', ChansForICA, 'pca',EffectiveRank);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2, 'setname', [SUB{i} '_ds_addChans_PREP_bp_refs_event_icaPrep2_weighted'], 'savenew', [Subject_Path SUB{i} '_ds_addChans_PREP_bp_refs_event_icaPrep2_weighted.set'], 'gui', 'off');
        
        %Load the continuous EEG data file outputted from Script #1 in .set EEGLAB file format
        EEG = pop_loadset( 'filename', [SUB{i} '_ds_addChans_PREP_bp_refs_event.set'], 'filepath', Subject_Path);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3, 'setname', [SUB{i} '_ds_addChans_PREP_bp_refs_event'], 'gui', 'off');
        
        %Transfer ICA weights to the continuous EEG data file (e.g., without the break periods and noisy segments of data removed)
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4, 'setname', [SUB{i} '_ds_addChans_PREP_bp_refs_event_icaWeighted'], 'savenew', [Subject_Path SUB{i} '_ds_addChans_PREP_bp_refs_event_icaWeighted.set'], 'gui', 'off');
        ALLEEG(CURRENTSET).icachansind = ALLEEG(2).icachansind;
        ALLEEG(CURRENTSET).icaweights = ALLEEG(2).icaweights;
        ALLEEG(CURRENTSET).icawinv = ALLEEG(2).icawinv;
        ALLEEG(CURRENTSET).icasphere = ALLEEG(2).icasphere;
        ALLEEG(CURRENTSET).chaninfo.icachansind = ALLEEG(2).chaninfo.icachansind;
        % make sure current EEG is populated with this new info
        EEG = ALLEEG(CURRENTSET);
        
        % new, an automated ICA removal using ICLabel plug in.
        % Using ICLabel output is more about removing components.
        % So it is placed in  process X4 instead.
        EEG = iclabel(EEG, 'default');
        
        % save the newly complete file to disk.
        EEG = pop_saveset( EEG, 'filename', [SUB{i} '_ds_addChans_PREP_bp_refs_event_icaWeighted.set'],'filepath', Subject_Path);

    end % End subject loop

end
%***********************************************************************************************************************************************
