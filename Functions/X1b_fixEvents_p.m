
% updated for multithread 17.08.21
% updated to allow relocating file to Functions folder. 18.08.21

% updated 15.09.21 to print event code counts to data directory, as some 
% files are crashing later because missing event codes. Best to 
% catch this early. Prints .csv for raw values and a separate .csv for
% trigger values that are > 256 and can be corrected to 0-255 range. 

function X1b_fixEvents_p(DataConfig, SUB)
% initialize
close all;
clearvars -except DataConfig SUB;

    % Location of the main study directory
    DIR = fileparts(pwd);
    
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
            
            if isempty(DataConfig.EventsToAdjust)
                % then don't adjust any events.
            else % else shift the specified events by the specified amount.
                
                % check to see that first nominated event was moved as
                % intended. (Assume the rest follow ok).
                preEvents = struct2table(EEG.event);
                % sometimes this will change edftype into cell of strings,
                % sometimes into cell of numerics. Need to force string.
                if strcmp(class(preEvents.edftype{1}), 'double')
                    % it's made a cell of numbers, so force a change.
                    preEvents.edftype = cellfun(@num2str, preEvents.edftype, 'UniformOutput', false);
                end
                indxtmp = strcmp(preEvents.type, DataConfig.EventsToAdjust{1});
                preEvents = preEvents(indxtmp,:);
                
                % now adjust the specified events.
                % should just use the built in EEGlab function.
                EEG = pop_adjustevents(EEG,'addms', DataConfig.TimeShift{1}*1000, 'eventtypes', DataConfig.EventsToAdjust);
                % now check the adjustment was sound. 
                postEvents = struct2table(EEG.event);
                % sometimes this will change edftype into cell of strings,
                % sometimes into cell of numerics. Need to force string.
                if strcmp(class(postEvents.edftype{1}), 'double')
                    % it's made a cell of numbers, so force a change.
                    postEvents.edftype = cellfun(@num2str, postEvents.edftype, 'UniformOutput', false);
                end
                indxtmp = strcmp(postEvents.type, DataConfig.EventsToAdjust{1});
                postEvents = postEvents(indxtmp,:);
                
                % compare the pre and post events and output the change.
                heightDiff = height(preEvents) - height(postEvents);
                if heightDiff == 0
                    % then they're the same length, so just compare
                    % directly.
                    latencyDiff = preEvents.latency - postEvents.latency;
                    % initially measured in samples, so convert to sec.
                    latencyDiff = latencyDiff./EEG.srate;
                    plot(latencyDiff);
                    title('AdjustmentPerEvent in seconds');
                    saveas(gcf, [Subject_Path filesep 'X1b_' SUB{i} '_EventsMoved.png']);
                else if heightDiff > 0
                        % then one or more events were chopped off pre to
                        % post.
                        preEvents = preEvents(heightDiff+1:end, :);
                    latencyDiff = preEvents.latency - postEvents.latency;
                    % initially measured in samples, so convert to sec.
                    latencyDiff = latencyDiff./EEG.srate;
                    plot(latencyDiff);
                    title('AdjustmentPerEvent in seconds');
                    saveas(gcf, [Subject_Path filesep 'X1b_' SUB{i} '_EventsMoved.png']);
                    else % something has gone wrong. Can't be adding events.
                        disp('Error: Somehow events were added when shifting times.');
                    end
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