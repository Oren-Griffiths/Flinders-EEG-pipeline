
% little function to force the configuration data to the right format, size
% etc needed by the sub-functions. Prior to this call, all fields of
% DataConfig are strings. After this call, most fields are cells that
% contain either scalars or strings. The last few entries to do with
% identifying particular channels are not cells; they are doubles. 

% last updated 18.08.21

function DataConfig = adjustConfigData(DataConfig)

try
    % adjust SUB. Is a char (or empty char). Needs to be an empty cell, a
    % 1x1 cell (with char inside) or a 1 x X cell (with chars inside).
    if ~isempty(DataConfig.SUB)
        if ~isempty(strfind(DataConfig.SUB, ';'))
            % more than one value.
            DataConfig.SUB = strsplit(DataConfig.SUB, ';');
            % remove any leading/trailing whitspace. Saves errors.
            for k = 1:length(DataConfig.SUB)
                DataConfig.SUB{k} = strtrim(DataConfig.SUB{k});
            end
        else % just one subject.
            % remove any leading/trailing whitspace. Saves errors.
            DataConfig.SUB = cellstr(strtrim(DataConfig.SUB));
        end
    else
        % won't trigger an error now, but will fail as soon as subfunction
        % is called a folder is tried to be found.
        DataConfig.SUB = {};
    end
    
    % BinListing. Take string, make cell of string.
    % same for RawFileType
    DataConfig.BinListing = cellstr(strtrim(DataConfig.BinListing));
    DataConfig.RawFileType = cellstr(strtrim(DataConfig.RawFileType));
    DataConfig.ReReference = cellstr(strtrim(DataConfig.ReReference));
    DataConfig.PREP = cellstr(strtrim(DataConfig.PREP));
    % TotalChannels. Take string, make numeric, place in cell.
    % same for DownSample
    % Same for BadChanThreshold
    % same for timeshift.
    % same for rawSrate
    % same for srate
    DataConfig.TotalChannels = num2cell(str2num(DataConfig.TotalChannels));
    DataConfig.DownSample = num2cell(str2num(DataConfig.DownSample));
    DataConfig.BadChanThreshold = num2cell(str2num(DataConfig.BadChanThreshold));
    DataConfig.TimeShift = num2cell(str2num(DataConfig.TimeShift));
    DataConfig.rawSrate = num2cell(str2num(DataConfig.rawSrate));
    DataConfig.HPfilter = num2cell(str2num(DataConfig.HPfilter));
    DataConfig.LPfilter = num2cell(str2num(DataConfig.LPfilter));
    DataConfig.FiltOrder = num2cell(str2num(DataConfig.FiltOrder));
    DataConfig.RemoveBlinks = num2cell(str2num(DataConfig.RemoveBlinks));
    DataConfig.CustomAR =  num2cell(str2num(DataConfig.CustomAR));
%     if ~isempty(DataConfig.srate)
%         DataConfig.srate = num2cell(str2num(DataConfig.srate)); % not a typo; current srate is rawSrate if loading the file new.
%     else % if it's empty put in the raw sampling rate
%         DataConfig.srate = DataConfig.rawSrate;
%     end
    DataConfig.EpochMin = num2cell(str2num(DataConfig.EpochMin));
    DataConfig.EpochMax = num2cell(str2num(DataConfig.EpochMax));
    DataConfig.BaselineMin = num2cell(str2num(DataConfig.BaselineMin));
    DataConfig.BaselineMax = num2cell(str2num(DataConfig.BaselineMax));
    
    % do a quick check and reorder mins and maxes if necessary
    if isempty(DataConfig.EpochMin) || isempty(DataConfig.EpochMin)
    else
        if DataConfig.EpochMax{1} < DataConfig.EpochMin{1}
            % do the old switcharoo to make the min into the max.
            tVar = DataConfig.EpochMax;
            DataConfig.EpochMax = DataConfig.EpochMin;
            DataConfig.EpochMin = tVar;
        end
    end
    
    % do a quick check and reorder mins and maxes if necessary
    if isempty(DataConfig.BaselineMin) || isempty(DataConfig.BaselineMin)
    else
        if DataConfig.BaselineMax{1} < DataConfig.BaselineMin{1}
            % do the old switcharoo to make the min into the max.
            tVar = DataConfig.BaselineMax;
            DataConfig.BaselineMax = DataConfig.BaselineMin;
            DataConfig.BaselineMin = tVar;
        end
    end
    
    % "events to adjust". Is a char (or empty char). Needs to be an empty cell, a
    % 1x1 cell (with char inside) or a 1 x X cell (with chars inside).
    if ~isempty(DataConfig.EventsToAdjust)
        % there's something there to ajust.
        if ~isempty(strfind(DataConfig.EventsToAdjust, ';'))
            % 2 or more properly formatted entries.
            DataConfig.EventsToAdjust = strsplit(DataConfig.EventsToAdjust, ';');
            for k = 1:length(DataConfig.EventsToAdjust)
                DataConfig.EventsToAdjust{k} = strtrim(DataConfig.EventsToAdjust{k});
            end
        else % 1 entry. Make it numeric and save it in a cell.
            DataConfig.EventsToAdjust = cellstr(strtrim(DataConfig.EventsToAdjust));
        end
    else % nothing there, just switch it to an empty cell.
        DataConfig.EventsToAdjust = {};
    end
    
%     % Bins. Is a char (or empty char). Needs to be an empty cell, a
%     % 1x1 cell (with num inside) or a 1 x X cell (with nums inside).
%     if ~isempty(DataConfig.Bins)
%         % there's something there to ajust.
%         if ~isempty(strfind(DataConfig.Bins, ';'))
%             % 2 or more properly formatted entries.
%             DataConfig.Bins = strsplit(DataConfig.Bins, ';');
%             for k = 1:length(DataConfig.Bins)
%                 DataConfig.Bins{k} = str2num(DataConfig.Bins{k});
%             end
%         else % 1 entry. Make it numeric and save it in a cell.
%             DataConfig.Bins = num2cell(str2num(DataConfig.Bins));
%         end
%     else % nothing there, just switch it to an empty cell.
%         DataConfig.Bins = {};
%     end
    
    % RelevantCodes" variable.
    % Is a char (or empty char). Needs to be an empty cell, a
    % 1x1 cell (with num inside) or a 1 x X cell (with nums inside).
    if ~isempty(DataConfig.RelevantCodes)
        % there's something there to adjust.
        if ~isempty(strfind(DataConfig.RelevantCodes, ';'))
            % 2 or more properly formatted entries.
            DataConfig.RelevantCodes = strsplit(DataConfig.RelevantCodes, ';');
            for k = 1:length(DataConfig.RelevantCodes)
                DataConfig.RelevantCodes{k} = str2num(DataConfig.RelevantCodes{k});
            end
        else % 1 entry. Make it numeric and save it in a cell.
            DataConfig.RelevantCodes = num2cell(str2num(DataConfig.RelevantCodes));
        end
    else % nothing there, just switch it to an empty cell.
        DataConfig.RelevantCodes = {};
    end
    
    % adjust CustomEpochs.
    % Is a char (or empty char). Needs to be an empty cell, a
    % 1x1 cell (with char inside) or a 1 x X cell (with chars inside).
    if ~isempty(DataConfig.CustomEpochs)
        if ~isempty(strfind(DataConfig.CustomEpochs, ';'))
            % more than one value.
            DataConfig.CustomEpochs = strsplit(DataConfig.CustomEpochs, ';');
            % remove any leading/trailing whitspace. Saves errors.
            for k = 1:length(DataConfig.CustomEpochs)
                DataConfig.CustomEpochs{k} = strtrim(DataConfig.CustomEpochs{k});
            end
            
        else
            DataConfig.CustomEpochs = cellstr(strtrim(DataConfig.CustomEpochs));
        end
    else
        % won't trigger an error now, but will fail as soon as subfunction
        % is called a folder is tried to be found.
        DataConfig.CustomEpochs = {};
    end
    
    %% Some parameters are baked in once you know the number of channels.
    % DataConfig.MaxBin{1} = max(cell2mat(DataConfig.Bins));
    
    % keep track of which process has been run, and what the name of the
    % last outputted file is... which will be nothing so far.
    
    DataConfig.LastSuffix = {};
    DataConfig.LastProcess = {};
    
    if DataConfig.TotalChannels{1} == 32
        DataConfig.KeyChans = num2cell([1 32 33 34 35 36 37 1 41]); % for 32 channel recording
        DataConfig.ChanSetup = cellstr('ChannelsFor32wMastoids.txt');
        DataConfig.AddCorrVEOG = cellstr('Add_ICAcorr_VEOGs_32.txt');
        DataConfig.ChanLocs = cellstr('standard-10-5-cap385.elp');
        DataConfig.rawVEOG = 39;
        DataConfig.rawHEOG = 38;
        DataConfig.corrVEOG = 43;
        DataConfig.corrHEOG = 42;
        DataConfig.cz_chan = 32;
        DataConfig.firstScalp = 1;
        DataConfig.lastScalp = 32;
        DataConfig.o1 = 15;
        DataConfig.o2 = 17;
        DataConfig.oz = 16;
        DataConfig.Fp1 = 1;
        % First scalp chan, Last scalp chan, L mastoid, R mastoid, HEOG L,
        % HEOG R, VEOG Low, VEOG high, lastChan]
    else % 64 channel recording
        DataConfig.ChanSetup = cellstr('ChannelsFor64wMastoids.txt');
        DataConfig.KeyChans = num2cell([1 64 65 66 67 68 69 1 73]);
        DataConfig.AddCorrVEOG = cellstr('Add_ICAcorr_VEOGs_64.txt');
        DataConfig.ChanLocs = cellstr('standard-10-5-cap385.elp');
        DataConfig.rawVEOG = 71;
        DataConfig.rawHEOG = 70;
        DataConfig.corrVEOG = 75;
        DataConfig.corrHEOG = 74;
        DataConfig.cz_chan = 51;
        DataConfig.firstScalp = 1;
        DataConfig.lastScalp = 64;
        DataConfig.o1 = 27;
        DataConfig.o2 = 64;
        DataConfig.oz = 29;
        DataConfig.Fp1 = 1;
    end
    
catch ME
    display('Error in config file. Workspace saved.');
    save('Debug_workspace.mat');
    rethrow(ME);
end

end