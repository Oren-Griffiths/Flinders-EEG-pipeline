%Script #5

%Operates on individual subject data
%Uses the output from Script #4: Remove_ICA_Components.m
%This script loads the semi-continuous ICA-corrected EEG data file from Script #4, creates an Event List containing a record of all event codes and their timing, assigns events to bins using Binlister, epochs the EEG, and performs baseline correction.

function X5_binEpochs_p(DataConfig, SUB)

% initialize
close all;
clearvars -except DataConfig SUB;

    % Location of the main study directory
    DIR = fileparts(pwd)
    
    % location of preprocessing files.
    Current_File_Path = pwd;
    
    %**********************************************************************************************************************************************************************
    
    %Loop through each subject listed in SUB
    for i = 1:length(SUB)
        
        %Open EEGLAB and ERPLAB Toolboxes
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        
        %Define subject path based on study directory and subject ID of current subject
        Subject_Path = [DIR filesep SUB{i} filesep];
        
        % Load the semi-continuous ICA-corrected EEG data file outputted from Script #4 in .set EEGLAB file format
        FileToOpen = [SUB{i} DataConfig.LastSuffix{1}];
        EEG = pop_loadset( 'filename', FileToOpen, 'filepath', Subject_Path);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [SUB{i} '_ds_PREP_ica_corr_cbip'], 'gui', 'off');
        
        
        if size(DataConfig.CustomEpochs) > 0 % custom text epochs to create.
            EventLabels = DataConfig.CustomEpochs;
            DataConfig.CustomEventTriggers = {}; % initialize output cell array.
            for ThisLabel = 1:numel(EventLabels)
                % create a token numeric label per text label. Start at 1001.
                DataConfig.CustomEventTriggers{ThisLabel} = 1000+ThisLabel;
                % cycle through labels and find their appearance
                for ThisEvent = 1:numel(EEG.event)
                    if strcmp (EEG.event(ThisEvent).type,EventLabels(ThisLabel))
                        EEG.event(ThisEvent).type = num2str(1000+ThisLabel);
                    end
                end %cycling vertically through candidate events.
            end % label by label loop
        end
        
        %Create EEG Event List containing a record of all event codes and their timing
        EEG  = pop_creabasiceventlist( EEG , 'AlphanumericCleaning', 'on', 'BoundaryNumeric', { -99 }, 'BoundaryString', { 'boundary' }, 'Eventlist', [Subject_Path SUB{i} '_Eventlist.txt'] );
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2, 'setname', [SUB{i} '_ds_PREP_ica_corr_cbip_elist'], 'savenew', [Subject_Path SUB{i} '_ds_PREP_ica_corr_cbip_elist.set'], 'gui', 'off');
        
        % must be using numeric codes with Binlister file.
        %Assign events to bins with Binlister; an individual trial may be assigned to more than one bin (bin assignments can be reviewed in each subject's '_Eventlist_Bins.txt' file)
        EEG  = pop_binlister( EEG , 'BDF', [Current_File_Path filesep 'SupportingDocs' filesep  DataConfig.BinListing{1}], 'ExportEL', [Subject_Path SUB{i} '_Eventlist_Bins.txt'], 'IndexEL',  1, 'SendEL2', 'EEG&Text', 'UpdateEEG', 'on', 'Voutput', 'EEG' );
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3, 'setname', [SUB{i} '_ds_PREP_ica_corr_cbip_elist_bins'], 'savenew', [Subject_Path SUB{i} '_ds_PREP_ica_corr_cbip_elist_bins.set'], 'gui', 'off');
        
        %Epoch the EEG into 1-second segments time-locked to the response (from -200 ms to 800 ms) and perform baseline correction using the average activity from -200 ms to 0 ms
        EEG = pop_epochbin( EEG , [DataConfig.EpochMin{1} DataConfig.EpochMax{1} ],  [ DataConfig.BaselineMin{1}, DataConfig.BaselineMax{1}]);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4, 'setname', [SUB{i} '_ds_PREP_ica_corr_cbip_elist_bins_epoch'], 'savenew', [Subject_Path SUB{i} '_ds_PREP_ica_corr_cbip_elist_bins_epoch.set'], 'gui', 'off');
        
        close all;
        
    end % End subject loop
    
%     DataConfig.LastProcess = cellstr('X5_BinEpochs');
%     DataConfig.LastSUB = SUB(i); % last participant processed.
%     DataConfig.LastSuffix = cellstr('_ds_PREP_ica_corr_cbip_elist_bins_epoch.set');
    
end % end of function.
%**********************************************************************************************************************************************************************
