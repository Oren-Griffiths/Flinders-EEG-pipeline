% global goal of this file is to use the same config file, add a few
% details to specify what kind of ERP you want, and then generates ERPs,
% figures, etc.

%% Get the extra details from the user [i.e. change these values.]
% what's the relevant config file called?
ConfigFileName = 'WIMR_Config_testing';

% 10:20 names of the channels you want. If you select more than one
% channel, it will average across them (i.e. treat it as a single montage).
% If you want to compare different channels/AOIs, run this script more
% than once with different channels chosen each time.
keyChans = {'FCz', 'Fz', 'Cz', 'FC1', 'FC2'};

% choose a time period you want to take the average across. Measured in ms.
measureWindow = [-1000, 1000];

% choose a frequency target (targetHz), and a visualization window 
% (rangeHz) around that value. 
targetHz = 15;
rangeHz = [3, 40];
% how do you want to do baselining?
baselining = 'none'; % or 'none' or 'z'

% if a person has less than X clean epochs are AR-rejection, then remove
% them from the averaging process. Set X. Put empty '[]' to ignore min.
minEpochs = 10;

% include the option of generating a "mask" file which has an 1/0 entry for
% every subject: 1s being included in the grand average analysis and 0s
% being excluded. If you just want to include everyone, write 'none'.
% default file created by X0 file is called 'PIDmask.xlsx'
maskFile = 'none';

% if you want to e.g. compare bin 1 with bin 3 in a 6 bin experiment, then
% put in a vector of:
% [1 0 -1 0 0 0 ]
% needs to have same number of entries as there are bins. Complex contrasts
% supported, e.g. [ 1 -0.5 -0.5 0 0 0]. If you just want all bins considered
% separately, leave it blank. Must be normalized (i.e. sum to 0), ...
% and ideally abs(sum) = 2 as well.
binContrast = [];

% do you want the output figures to show information about peak value and
% latency (in ERP waveforms) and min/max channel values in the topoplots?
% if you want this info, set variable to 1. Else leave as 0.
showPeakInfo = 1;
% do you want to add a fine overlay of each individual to the grand average
% waveforms (e.g. to check for presence of outliers)? If so, set to 1.
showIndividualTraces = 1;
% Ok, what does it do with this info?
% Generates a global average figure, a data set with raw sample-by-sample
% data per person, averaged values per participant during the measurement
% window. Places output in current folder. File is 'ERP_output.mat'

%% open the config file to grab rest of the relevant info.
% let's loop through all relevant files and draw the ERG1 figures.
Current_File_Path = pwd;
addpath('Functions');
ConfigFilePath = [Current_File_Path filesep 'SupportingDocs' filesep ConfigFileName '.xlsx'];
Options = detectImportOptions(ConfigFilePath);

for k = 1:numel(Options.VariableTypes)
    Options.VariableTypes{k} = 'char';
end
DataConfig = table2struct(readtable(ConfigFilePath, Options));
DataConfig = adjustConfigData(DataConfig);

% and open eeglab to access the EEGlab functions
eeglab;

%% ok, need to figure out some key values for data aggregation.

% find the number of bins direct from the binlister document.
filename = [pwd filesep 'SupportingDocs' filesep DataConfig.BinListing{1}];
fileID = fopen(filename);
BinListText = fscanf(fileID,'%s');
fclose(fileID);
temp = strfind(BinListText,'Bin');
NoOfBins = length(temp); % value we need.

% quickly open the first file to get the number of chans.
testFile_mat = [DataConfig.SUB{1} '_ARcorrectedBins.mat'];
testFolder= [fileparts(pwd) filesep DataConfig.SUB{1}];
load([testFolder filesep testFile_mat]);
for ThisBin = 1:numel(GoodTrials)
    if isempty(GoodTrials(ThisBin).data)
    else
        arraySizes = size(GoodTrials(ThisBin).data);
        NoOfChans = arraySizes(1); % value we need.
        LengthOfEpoch = arraySizes(2); % value we need.
        if length(arraySizes) > 2
            NoOfEpochs = arraySizes(3); % value we need.
        else
            NoOfEpochs = 1;
        end
        ThisBin = numel(GoodTrials); % enter the exit value.
    end
    
    % have you already grabbed the chanloc info?
    if exist('chanlocs')
        % already got the chanloc info. Ignore.
    else
        if ~isempty(GoodTrials(ThisBin).chanlocs)
            chanlocs = GoodTrials(ThisBin).chanlocs;  % we need this structure too.
        end
    end
end

% and let's get the values from another method too. Makes sure everything
% aligns.
LengthOfEpoch_2 = (DataConfig.EpochMax{1} - DataConfig.EpochMin{1})/1000 * ...
    DataConfig.DownSample{1};
srate = DataConfig.DownSample{1}; % need this value too.
if LengthOfEpoch_2 == LengthOfEpoch
    % nice match. Phew!
else % possible error.
    display('Possible error in epoch lengths.');
    display(['DataConfig says epoch should be this long:' num2str(LengthOfEpoch_2)]);
    display(['Outputted data are this long:' num2str(LengthOfEpoch)]);
end

% check that the values inputted are sensible. 
if targetHz < rangeHz(1) || targetHz > rangeHz(2)
    disp('****************************************************');
    disp('Target frequency value is not in the frequency range. ');
    disp('****************************************************');
end

% generate an x-axis for plotting. Measurement is now in seconds.
times = ([1:LengthOfEpoch].*1/srate) + (DataConfig.EpochMin{1}/1000);
% review the period that we wish to evaluate. 
measurePeriod = (times > measureWindow (1)/1000 & times < measureWindow(2)/1000);
times = times(measurePeriod);
LengthOfEpoch = length(times); % need to trim time series pre FFT (only in FFT).
FFT_length = floor(LengthOfEpoch/2)+1;
HzBins = 0: (srate/2)/(FFT_length-1)  : srate/2;

% initialize a variable full of not-numbers (0s,1s would bias means if an
% error was made somewhere down the line. Here errors make it break).
% structure: a n-element cell array (cells are bins). Each bin contains a
% 3d structure: PIDs by chans by samples.
for ThisBin = 1:NoOfBins
    participantAverages{ThisBin} = NaN(length(DataConfig.SUB),NoOfChans,FFT_length); % 3d data structure per bin.
end

% find out who to include in analysis and who to exclude.
if strcmp(maskFile, 'none')
    % no mask vector requested
    disp('No mask vector requested.');
    maskVector = [];
else
    % find the vector of INCs/EXCs
    if exist(['SupportingDocs' filesep maskFile])
        % read in the SUB IDs and the 1/0 status
        T = readtable(['SupportingDocs' filesep maskFile]);
        % convert that into a vector.
        if height(T) == length(DataConfig.SUB)
            disp('Mask vector requested, found and applied.');
            maskVector = T{:,2};
        else
            disp('Error in mask file. Mistmatched lengths. Including everyone.');
            maskVector = [];
        end
        
    else % no file available, so no mask vector generated.
        disp('No mask vector file found. None applied. All included.');
        maskVector = [];
    end
end


%% loop through SUBS and gather *per participant* averages by averaging across epochs (within bins).
SUB = DataConfig.SUB;
for k = 1:length(SUB)
    
    % just working space to check how data is binned in X7 code.
    testFolder= [fileparts(pwd) filesep SUB{k}];
    testFile_mat = [SUB{k} '_ARcorrectedBins.mat'];
    
    % load the outputted data in matlab format.
    load([testFolder filesep testFile_mat]);
    
    % adjust the data according to bin contrasts to streamline the data
    % down to a single (e.g. difference) time series, if that is what's
    % declared.
    
    if isempty(binContrast)
        % do nothing. leave the input alone.
    else % possible some bins are missing in some participants. correct for that.
        residual = length(binContrast) - numel(GoodTrials);
        % fill GoodTrials variable in with blank cells if needed.
        if residual == 0
            if residual < 0
                display('More bins in data than in binContrast vector. How?');
            else % add some empty bins into the GoodTrials variable.
                for missingBin = 1:residual
                    GoodTrials(missingBin+numel(GoodTrials)).data = [];
                    GoodTrials(missingBin+numel(GoodTrials)).ID = [];
                    GoodTrials(missingBin+numel(GoodTrials)).chanlocs = GoodTrials(1).chanlocs;
                    GoodTrials(missingBin+numel(GoodTrials)).srate = GoodTrials(1).srate;
                end
            end
            % the math/adjustment is done down below.
        end
    end
    
    
    for ThisBin = 1:NoOfBins
        if ThisBin > numel(GoodTrials)
            % do nothing, this entry is already all NaNs
        else % check to see if there are any data there.
            if isempty(GoodTrials(ThisBin).data)
                % do nothing, this entry is already all NaNs
                display(['Skipping SUB ' SUB{k} ' Bin ' num2str(ThisBin) '. No data.']);
            else
                if isempty(minEpochs) || size(GoodTrials(ThisBin).data,3) > minEpochs
                    
                    % just use the scalp channels for this analysis.
                    for ThisChan = 1:NoOfChans
                        % epochs by "samples" (output of FFT, so spectral
                        % powers).
                        disp(['FFT for Chan: ' num2str(ThisChan) ' for Bin:' num2str(ThisBin) ' PID: ' SUB{k}]);
                        % intialize output variable.
                        NoOfEpochs = size(GoodTrials(ThisBin).data,3);
                        FFT_placeholder = NaN(NoOfEpochs,FFT_length);
                        for ThisEpoch = 1:NoOfEpochs
                            switch baselining
                                case 'subtract'
                                    FFT_placeholder(ThisEpoch,:) = ...
                                        applyFFTbaseline(SimpleFFT(GoodTrials(ThisBin).data(ThisChan,measurePeriod,ThisEpoch)))';
                                case 'z'
                                    FFT_placeholder(ThisEpoch,:) = ...
                                        applyFFTbaselineZ(SimpleFFT(GoodTrials(ThisBin).data(ThisChan,measurePeriod,ThisEpoch)))';
                                case 'none'
                                    FFT_placeholder(ThisEpoch,:) = ...
                                        SimpleFFT(GoodTrials(ThisBin).data(ThisChan,measurePeriod,ThisEpoch));
                            end
                        end % of epoch by epoch loop.
                        % structure of data per bin: PIDs by chans by samples.
                        participantAverages{ThisBin}(k,ThisChan,:) = nanmean(FFT_placeholder,1);
                    end % channel by channel loop.
                end % of checking against min number of epochs. 
            end % of skipping empty data sets 2
        end % of skipping empty data sets 1
    end % of bin by bin loop.
end % of subject by subject loop

% sanity check to make sure it's outputting approriate spectral values.
% plot(HzBins,nanmean(FFT_placeholder));

%% correct or adjust the data, if necessary
% will do the bin contrast below when dimensionality reduced.
if ~isempty(binContrast)
    involvedBins = find(binContrast);
    for i = 1:length(involvedBins)
        if i == 1
            contrastAverages = participantAverages{involvedBins(i)}(:,:,:) .*binContrast(involvedBins(i));
        else
            contrastAverages = contrastAverages + participantAverages{involvedBins(i)}(:,:,:) .*binContrast(involvedBins(i));
        end
    end
    % only one bin now.
    NoOfBins = 1;
    % move the data back into four dimensional structure.
    clear participantAverages; % get rid of structure of variable. restate.
    participantAverages{1} = contrastAverages;
    clear contrastAverages; % big variable. ditch it when done.
end

% apply the masking vector to exclude some participants ,if appropriate.
if isempty(maskVector)
    % do nothing and leave participantAverages intact.
else
    for ThisBin = 1:NoOfBins
        participantAverages{ThisBin}(~maskVector, :, :) = NaN;
        % using NaNs keeps the number of values in that dimension the same, and
        % all aggregate statistics (mean, SD) exclude those values anyway.
    end
end

%% and now start to calculate the output the data needed.
outputFolder = 'FFT_GrandAverages';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

keyFreqs = (HzBins > rangeHz(1) & HzBins < rangeHz(2));
% a=[34.8 31.2 29 26.7 39.5];%dummy data
% n=33;
distFromTarget = abs(HzBins-targetHz);
[minDev,targetHz_idx] = min(distFromTarget);
targetHz_actual = HzBins(targetHz_idx);

% get channel numbers
% that is, convert channel indices from channel names.
% have done this assuming that all bins have the same channel location
% structure, but robust to missing values for some bins.

for ThisChan = 1:length(keyChans)
    keyChanIdx(ThisChan) = find(strcmp({chanlocs.labels}, keyChans{ThisChan})==1);
    ThisBin = NoOfBins; % skip out of the loop.
end

% participantAverages.
% structure: a n-element cell array (cells are bins). Each bin contains a
% 3d structure: PIDs by chans by samples.
% average across channels (if more than one key channel chosen).
for ThisBin = 1:NoOfBins
    % pool the channels (if you have a multi channel montage).
    temp = squeeze(nanmean(participantAverages{ThisBin}(:,keyChanIdx,:),2));
    % now in PID by samples (freqs) format.
    tempForOutput{ThisBin} = temp;
end

% save that raw data (i.e. tempForOutput)
rawdataFilename = [pwd filesep outputFolder filesep 'RawOutput_PIDbyHz.xlsx'];
for ThisBin = 1:NoOfBins
    % declare a useful tab name
    if isempty(binContrast)
        tabname = ['Bin' num2str(ThisBin)];
    else
        tabname = 'DifferenceFFT';
    end
    % write the data
    % limit to the relevant epoch Hz range.
    writematrix(tempForOutput{ThisBin}(:,keyFreqs), rawdataFilename, 'Sheet', tabname, 'Range','B3');
    % write the row headers
    writecell(DataConfig.SUB', rawdataFilename, 'Sheet', tabname, 'Range', 'A3');
    writecell({'Hz'}, rawdataFilename, 'Sheet', tabname, 'Range', 'A1');
    writecell({'PID'}, rawdataFilename, 'Sheet', tabname, 'Range', 'A2');
    % write the column headers
    writematrix(HzBins(keyFreqs), rawdataFilename, 'Sheet', tabname, 'Range','B1');
    % and no need do have specific sheet for key DV, as it's contained
    % here. Just use the *** indicator to mark closest value to target
    % frequency so it's easy to tell.
    emptyCell = cell(1,length(HzBins)); 
    emptyCell{targetHz_idx} = '***';
    writecell(emptyCell(keyFreqs), rawdataFilename, 'Sheet', tabname, 'Range','B2');
end


%% and now start drawing.
% participantAverages.
% structure: participants, chans, samples, by bins.

% we want global mean and SEM over the keyPeriod time window.
% and we want separate figures for each bin.

for ThisBin = 1:NoOfBins
    HzBinsToPlot = HzBins(keyFreqs);
    % report grand average(across PIDs) for nominated chans in this
    % bin.
    % data in PID by samples format, so average across PIDs.
    FreqsToPlot = nanmean(tempForOutput{ThisBin}(:,keyFreqs),1);
    
    % find peak within measurement window?
    % PIDs by chans by samples.
    [posPeak, posPeakIdx] = max(FreqsToPlot);
    posPeakHz = HzBinsToPlot(posPeakIdx);
    posPeakText = ['Max Val = ' num2str(round(posPeak,2)) 'uV^2/Hz at ' ...
        num2str(round(posPeakHz,2)) 'Hz' ];

    % start drawing
    figure;
    hold on
    if showIndividualTraces == 1
        for ThisPID = 1:size(tempForOutput{ThisBin},1)
            line(HzBinsToPlot, tempForOutput{ThisBin}(ThisPID,keyFreqs), 'LineStyle', ':', 'Color', 'k', 'LineWidth', 0.5);
        end
    end
    % main line to draw.
    line(HzBinsToPlot, FreqsToPlot, 'LineStyle', '-', 'Color', 'k', 'LineWidth', 2);
    hold off
    
    % extra details for figure.
    if isempty(binContrast)
        title(['FFT of target range in bin: ' num2str(ThisBin) ' at channel(s): ' string(strjoin(keyChans))]);
    else
        title(['FFT of target range in difference wave at channel(s): ' string(strjoin(keyChans))]);
    end
    if showPeakInfo == 1
        % provide some info overlays about peaks.
        text(posPeakHz, posPeak, ['\leftarrow ' posPeakText], 'Color','green','FontSize',10);
    end
    ylabel('Spectral power density(uV^2/Hz)');
    xlabel('Hz');
    y_cap = 2*max(abs(FreqsToPlot));
    ylim([-1*y_cap, y_cap]);
    % change some values and save.
    f = gcf;
    f.Units = 'inches';
    f.OuterPosition = [0.5 0.5 5.5 5.5]; % make the figure 5 inches in size.
    fig_filename = [outputFolder filesep 'Bin_' num2str(ThisBin) '_GrandAverageFFT.png'];
    disp(['Saving ERP image ' fig_filename]);
    exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
    close(gcf);
    
    % and now do a topoplot per bin.
    
    % here we just need to pull the value per channel associated with the
    % target frequency. We already know it's here: targetHz_idx
    % And the raw data to draw from is: participantAverages.
    % structure: a n-element cell array (cells are bins). Each bin contains a
    % 3d structure: PIDs by chans by samples.
    % average across channels (if more than one key channel chosen).
    
    dataToPlot = nanmean(participantAverages{ThisBin} (:,:,targetHz_idx), 1);
    % limit the plotting to scalp channels.
    dataToPlot = dataToPlot(1:DataConfig.TotalChannels{1});
    
    if sum(isnan(dataToPlot)) == length(dataToPlot)
        % don't actually draw the image if there aren't any data.
    else
        % and grab some information about min/max
        [minChanVal, minChanIdx] = min(dataToPlot);
        [maxChanVal, maxChanIdx] = max(dataToPlot);
        minChanLbl = chanlocs(minChanIdx).labels;
        maxChanLbl = chanlocs(maxChanIdx).labels;
        minText = ['Min Chan: ' minChanLbl ' at ' num2str(round(minChanVal,2)) 'uV'];
        maxText = ['Max Chan: ' maxChanLbl ' at ' num2str(round(maxChanVal,2)) 'uV'];
        figure;
        if showPeakInfo == 1
            % you've requested info about electrode locations.
            topoplot(dataToPlot, chanlocs, 'electrodes' ,'ptslabels');
        else
            % you just want a standard image. No info overlays.
            topoplot(dataToPlot, chanlocs);
        end
        colorbar;
        if isempty(binContrast)
            titleText = ['Topoplot of Bin: ' num2str(ThisBin) ...
                ' at ' num2str(targetHz_actual) 'Hz'];
        else
                        titleText = ['Topoplot of difference wave at ' ...
                num2str(targetHz_actual) 'Hz'];
        end
        title(titleText);
        colourCap = 1.5*max(abs(dataToPlot));
        caxis([-1*colourCap, colourCap]);
        if showPeakInfo == 1
            % add text related to min/max
            annotation('textbox', [0, 0.15, 0, 0], 'string', minText, 'FitBoxToText','on');
            annotation('textbox', [0, 0.05, 0, 0], 'string', maxText, 'FitBoxToText','on');
        end
        f = gcf;
        f.Units = 'inches';
        f.OuterPosition = [0.5 0.5 5.5 5.5]; % make the figure 5 inches in size.
        fig_filename = [outputFolder filesep 'Bin_' num2str(ThisBin) '_TargetHzTopoplot.png'];
        disp(['Saving topoimage image ' fig_filename]);
        exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
        close(gcf);
    end
end



