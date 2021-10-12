function X1b_fixEvents_p(DataConfig, SUB)
% updated for multithread 17.08.21
% updated to allow relocating file to Functions folder. 18.08.21

% updated 15.09.21 to print event code counts to data directory, as some 
% files are crashing later because missing event codes. Best to 
% catch this early. Prints .csv for raw values and a separate .csv for
% trigger values that are > 256 and can be corrected to 0-255 range. 


% initialize
close all;
clearvars -except DataConfig SUB;

    % Location of the main study directory
    DIR = fileparts(pwd)
    
    % location of preprocessing files.
    Current_File_Path = pwd;
    
    %% events to correct and how to adjust them.
    % some events need to be moved.
    
    % Key input and preprocess loop
    %Loop through each subject listed in SUB
    for i = 1:length(SUB)
        
        if isempty(DataConfig.EventsToAdjust)
            % no events to adjust so will be empty. So just update process
            % at very end, but don't open or save a file.
            
        else % some events to adjust.
            
            % Open EEGLAB and ERPLAB Toolboxes
            % Acts as initialization of the relevant variables per participant.
            [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
            
            % Define subject path based on study directory and subject ID of current subject
            Subject_Path = [DIR filesep num2str(SUB{i}) filesep];
            
            % record where we're up to in case of crash.
            DataConfig.CurrentSUB = SUB{i};
            
            % open the relevant file and get cracking.
            % should be an earlier file to load, but just in case there's not.
            if isempty(DataConfig.LastProcess)
                % this is the first process for this file, so open the raw data.
                switch DataConfig.RawFileType{1}
                    case '.bdf'
                        FileToOpen = [SUB{i} '.bdf'];
                        EEG = pop_biosig([Subject_Path  FileToOpen], 'ref', DataConfig.KeyChans{5});
                        EEG = eeg_checkset( EEG );
                    case '.xdf'
                        FileToOpen = [SUB{i} '.xdf'];
                        EEG = pop_loadxdf([Subject_Path  FileToOpen]);
                        EEG = eeg_checkset( EEG );
                end
            else
                % not the first file to open. So grab the last file that we
                % processed.
                FileToOpen = [SUB{i} DataConfig.LastSuffix{1}];
                EEG = pop_loadset( 'filename', FileToOpen, 'filepath', Subject_Path);
            end
            
            %Load the continuous EEG data file outputted from Script #1a in .set EEGLAB file format
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [SUB{i} '_ds_addChans_PREP_bp_refs'], 'gui', 'off');
            
            % So need to do something like this combo:
            for ThisEventType_idx = 1:numel(DataConfig.EventsToAdjust)
                ThisEventType = DataConfig.EventsToAdjust{ThisEventType_idx}; % a character indicating the relevant event label (often numeric label).
                event_index = find(strcmp({EEG.event.type}, ThisEventType) == 1); % find instances of that label.
                for InstanceOfEvent = 1:length(event_index) % loop through instances of that label.
                    EEG.event(InstanceOfEvent).latency = EEG.event(InstanceOfEvent).latency + ...
                        EEG.srate *DataConfig.TimeShift{1} * 1000; % adjust latency of that event.
                end % instance by instance loop
                DataConfig.EventsAdjusted{ThisEventType_idx} = sum(event_index); % record how many events were changed, ordered by type.
            end % event by event loop
            
            % remove any events that have a negative latency in case that interferes later.
            % just loop  through all of them as finding negatives doesn't
            % seem to work (in my hands).
            % initialize an output field so that you can capture how many
            % events were removed.
            DataConfig.RemovedNegativeEvents = num2cell(0);
            for ThisEvent = numel(EEG.event):-1:1 % reverse order loop.
                if EEG.event(ThisEvent).latency < 0
                    EEG.event(ThisEvent) = []; % delete that entry
                    DataConfig.RemovedNegativeEvents{1} = DataConfig.RemovedNegativeEvents{1} + 1;
                end
            end
            
            % format for creating a new .set file and saving it.
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2, 'setname', [SUB{i} '_ds_addChans_PREP_bp_refs_event'], 'savenew', [Subject_Path SUB{i} '_ds_addChans_PREP_bp_refs_event.set'], 'gui', 'off');
            
        end % of the part to skip if no events to adjust
        
        % and even if it is skipped, output a table of the number of events
        % in a .txt file. But exclude 0,1, 256 codes which are just noise.
        disp('Calculating and reporting event tallies.');
        eventTable = struct2table(EEG.event);
        % convert to numeric and add a column called "codes"
        eventCodes = str2double(eventTable.type);
        eventTable.codes = eventCodes;
        % filter 256, 0, 1, values. 
        nanCodes = isnan(eventCodes); % remove text codes
        badCodes = eventCodes == 1 | eventCodes ==0 | eventCodes ==256;
        % removes pointless 0, 1, 256 auto-generated codes.
        goodCodes = ~(nanCodes | badCodes);
        eventTable = eventTable(goodCodes,:);
        % remove unnecessary fields.
        eventTable = removevars(eventTable, {'type', 'edftype', 'latency', ...
             'urevent'});
        for ThisRow = 1:height(eventTable)
            if eventTable.codes(ThisRow) > 256 
                eventTable.corrCodes(ThisRow) = eventTable.codes(ThisRow) - 256;
            end
        end
        
        % count the unique values and report them as .csv files.
        [C,ia,ic] = unique(eventTable.codes);
        code_counts = accumarray(ic,1);
        countsPerCode = array2table([C, code_counts]);
        countsPerCode.Properties.VariableNames(1:2) = {'event','count'};
        outFileName = [Subject_Path SUB{i} '_eventCounts.csv'];
        writetable(countsPerCode, outFileName);
        % do the same, but with 256 adjustment. 
        [C,ia,ic] = unique(eventTable.corrCodes);
        code_counts = accumarray(ic,1);
        countsPerCorrCode = array2table([C, code_counts]);
        countsPerCorrCode.Properties.VariableNames(1:2) = {'event','count'};
        outFileName= [Subject_Path SUB{i} '_eventCountsWith256Correction.csv'];
        writetable(countsPerCorrCode, outFileName);
        
    end % of subject by subject loop
    
    % record the last process performed. (pointless now with parallel
    % execution as DataConfig is not a global or passed out).

end % of function.