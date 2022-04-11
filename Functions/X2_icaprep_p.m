%Script #2
%Operates on individual subject data
%Uses the output from Script #1: Import_Raw_EEG_Reref_DS_Hpfilt.m
%This script loads the outputted continuous EEG data file from Script #1, removes segments of EEG during the break periods in between trial blocks, and
%removes especially noisy segments of EEG during the trial blocks to prepare the data for ICA. Note that the goal of this stage of processing is to remove
%particularly noisy segments of data; a more thorough rejection of artifacts will be performed later on the epoched data.

% updated 17.08.21 to make compatible with multithreading
% updated 20.08.21 to gracefully manage datasets with zero of the nominated
% trigger codes. 

function X2_icaprep_p(DataConfig, SUB)

% initialize
close all;
clearvars -except DataConfig SUB;


    % Location of the main study directory
    DIR = fileparts(pwd);
    
    % location of preprocessing files.
    Current_File_Path = pwd;
    
    
    %*************************************************************************************************************************************
    % Loop through each subject listed in SUB
    for i = 1:length(SUB)
        
        %Open EEGLAB and ERPLAB Toolboxes
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        
        %Define subject path based on study directory and subject ID of current subject
        Subject_Path = [DIR filesep SUB{i} filesep];
        
        % record where we're up to in case of crash.
        DataConfig.CurrentSUB = SUB{i};
        
        % Load the continuous EEG data file outputted from Script #1 in .set EEGLAB file format
        FileToOpen = [SUB{i} DataConfig.LastSuffix{1}];
        EEG = pop_loadset( 'filename', FileToOpen, 'filepath', Subject_Path);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [SUB{i} '_ds_addChans_PREP_bp_refs_event'], 'gui', 'off');
        czPreClean = EEG.data(DataConfig.cz_chan, :); % grab dirty data for later plotting.
        
        %Remove segments of EEG during the break periods in between trial blocks (defined as 5 seconds or longer in between successive stimulus event codes)
        if ~isempty(DataConfig.RelevantCodes)
            % ok, but if there are no events consistent with the relevant
            % codes this will crash. So need to check for that too. 
            matchedEvents = 0; % start at 0 and update if a match is found. 
            % remove the extra fields (if present).
            tempStruct = EEG.event;
            if isfield(tempStruct, 'urevent')
                tempStruct = rmfield(tempStruct, 'urevent');
            end
            if isfield(tempStruct, 'latency')
                tempStruct = rmfield(tempStruct, 'latency');
            end
            if isfield(tempStruct, 'duration')
                tempStruct = rmfield(tempStruct, 'duration');
            end
            eventCodes = struct2cell(tempStruct);
            
            for k = 1:length(DataConfig.RelevantCodes) % loop through codes to check.
                if isempty(find( [eventCodes{:}] == DataConfig.RelevantCodes{k}  ))
                     % no matches found for this code. Do nothing. 
                else % bingo, found a match. Can do specific trimming. 
                    matchedEvents = 1; 
                end
            end
            
            if matchedEvents == 1
                EEG  = pop_erplabDeleteTimeSegments( EEG , 'displayEEG', 0, 'endEventcodeBufferMS',  500, 'ignoreUseEventcodes', cell2mat(DataConfig.RelevantCodes), 'ignoreUseType', 'Use', 'startEventcodeBufferMS',  500, 'timeThresholdMS',  2000 );
            else % else just general trimming.
                disp('ALERT: This dataset has 0 instances of the relevant triggers.');
                EEG  = pop_erplabDeleteTimeSegments( EEG , 'timeThresholdMS',  2000 );
            end
        else
            EEG  = pop_erplabDeleteTimeSegments( EEG , 'timeThresholdMS',  2000 );
        end
        
        %Load parameters for rejecting especially noisy segments of EEG during trial blocks from Excel file ICA_Prep_Values.xls. Default parameters can be used initially but may need
        % to be modified for a given participant on the basis of visual inspection of the data.
        
        % Excel file currently has sensible defaults (800, 2000, 100).
        [ndata, ~, alldata] = xlsread([Current_File_Path filesep 'SupportingDocs' filesep 'ICA_Prep_Values']);
        for j = 1:size(alldata,1)
            if isequal(SUB{i},num2str(alldata{j,1}))
                AmpthValue = alldata{j,2};
                WindowValue = alldata{j,3};
                StepValue = alldata{j,4};
            end
        end
        
        % Delete segments of the EEG exceeding the thresholds defined above
        EEG = pop_continuousartdet( EEG, 'ampth', AmpthValue, 'winms', WindowValue, 'stepms', StepValue, 'chanArray', DataConfig.KeyChans{1}:DataConfig.KeyChans{2}, 'review', 'off');
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3, 'setname', [SUB{i} '_ds_addChans_PREP_bp_refs_event_icaPrep2'], 'savenew', [Subject_Path SUB{i} '_ds_addChans_PREP_bp_refs_event_icaPrep2.set'], 'gui', 'off');
        czPostClean = EEG.data(DataConfig.cz_chan, :); % grab clean data for later plotting.
        
        % plot and save the degree of cleaning (stretching to match lengths)
        times1 = 1 : length(czPreClean);
        times1 = times1/EEG.srate;
        times2 = 1 : length(czPostClean);
        times2 = times2/EEG.srate;
        stretchfactor =  length(czPreClean) / length(czPostClean) ;
        plot(times1, czPreClean, stretchfactor*times2, czPostClean);
        txt = ['Amp = ', num2str(AmpthValue), ', Window = ', num2str(WindowValue), ', Step = ', num2str(StepValue)];
        ylim([-1000, 1000]);
        x_insert = max(times1)/2;
        y_insert = 1050;
        text(x_insert, y_insert, txt);
        save2pdf([Subject_Path filesep 'Figures' filesep 'X2_' SUB{i}  '_CleanedPreICA.pdf']);
        close all
        
    end  %End of subject loop

end
%*************************************************************************************************************************************
