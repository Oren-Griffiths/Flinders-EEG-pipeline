% global goal of this file is to use the same config file, add a few
% details to specify what kind of ERP you want, and then generates ERPs,
% figures, etc.

%% Get the extra details from the user [i.e. change these values.]
% what's the relevant config file called?
ConfigFileName = 'Config_Danielle_051022';

% 10:20 names of the channels you want. If you select more than one
% channel, it will average across them (i.e. treat it as a single montage).
% If you want to compare different channels/AOIs, run this script more
% than once with different channels chosen each time.
% for fronto-central montage
keyChans = { 'Cz', 'C3', 'C4',  'FC1', 'FC2', 'Fz', 'F3', 'F4' };
% for parietal sites (p3a)
% keyChans = { 'PO3', 'PO4', 'Pz',  'P3', 'P4', 'CP1', 'CP2' };

% what time period (ms) do you want visualized. This will obviously break
% if you choose an area large than the epoch declared in DataConfig, so
% choose sensibly.
wholeEpoch = [-200, 400];

% choose a time period you want to take the average across. Measured in ms.
measureWindow = [100, 200];

% if a person has less than X clean epochs are AR-rejection, then remove
% them from the averaging process. Set X. Put empty '[]' to ignore min.
minEpochs = 1;

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

% You can specify an extra lowpass filter that may be needed for terrible
% data. By default this should be left as 'none', but for messy data
% you may be forced to filter out HF noise again. You can do the same as
% what you've done previously (type 'same'). Otherwise choose a value e.g.
% 40Hz and enter it as a string '40'.
extraLowPass = 'none';

% add in a vector of how many samples to shift one or more bins (if, for
% example, there is a triggering issue where one trigger waits a frame or
% two and others do not). Should usually be empty. If there are values,
% then you'll need to supply a number of samples per bin (in a vector),
% e.g. for a 4 bin experiment where you want to shift bin 4 only [0,0,0,3];
sampleCorrection = [0,0,0,0,0,0,0,0,0,0,0,0];

% do you want the output figures to show information about peak value and
% latency (in ERP waveforms) and min/max channel values in the topoplots?
% if you want this info, set variable to 1. Else leave as 0.
showPeakInfo = 0;
% do you want to add a fine overlay of each individual to the grand average
% waveforms (e.g. to check for presence of outliers)? If so, set to 1.
showIndividualTraces = 1;
% Ok, what does it do with this info?
% Generates a global average figure, a data set with raw sample-by-sample
% data per person, averaged values per participant during the measurement
% window. Places output in current folder. File is 'ERP_output.mat'

% For the visualizations, do you want to have a fixed y-axis? 
% If not, leave it empty. 
fixedYlim = [-5,5];

% can declare the colour scheme in advance (for topos). 
colScheme = 'jet';

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

% get channel numbers
% that is, convert channel indices from channel names.
% have done this assuming that all bins have the same channel location
% structure, but robust to missing values for some bins.
for ThisChan = 1:length(keyChans)
    keyChanIdx(ThisChan) = find(strcmp({chanlocs.labels}, keyChans{ThisChan})==1);
    ThisBin = NoOfBins; % skip out of the loop.
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

% initialize a variable full of not-numbers (0s,1s would bias means if an
% error was made somewhere down the line. Here errors make it break).
% structure: a n-element cell array (cells are bins). Each bin contains a
% 3d structure: PIDs by chans by samples.
for ThisBin = 1:NoOfBins
    participantAverages{ThisBin} = NaN(length(DataConfig.SUB),NoOfChans,LengthOfEpoch); % 3d data structure per bin.
end

% generate an x-axis for plotting. Measurement is now in seconds.
times = ([1:LengthOfEpoch].*1/srate) + (DataConfig.EpochMin{1}/1000);

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

% figure out the extra lowpass filter, if required
if strcmp(class(extraLowPass), 'char')
    switch extraLowPass
        case 'none'
            extraLowPass_num = [];
            disp('No extra filter applied');
        case 'same'
            extraLowPass_num = DataConfig.LPfilter{1};
            disp(['Extra filter applied: ' num2str(DataConfig.LPfilter{1}) 'Hz' ]);
        otherwise
            extraLowPass_num = str2num(extraLowPass);
            disp(['Extra filter applied: ' extraLowPass 'Hz']);
    end
else
    disp('Error in supplied extra filter. None applied');
    extraLowPass_num = [];
end
% extraLowPass = DataConfig.LPfilter{1}; %


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
            else
                % check to see if there are enough epochs for a stable
                % estimate.
                if isempty(minEpochs) || size(GoodTrials(ThisBin).data,3) > minEpochs
                    % load the mean per epoch into a global variable.
                    temp = squeeze(nanmean(GoodTrials(ThisBin).data,3));
                    for ThisChan = 1:NoOfChans
                        if ~isempty(extraLowPass_num)
                            % particularly noisy data may need an extra
                            % lowpass filter (if cleanline + BP doesn't
                            % grab everything)
                            if k == 1 && ThisBin == 1
                                % you only need to design the filter on the
                                % first loop.
                                [~,dFilter] = lowpass(temp(ThisChan,:),extraLowPass_num, srate);
                                % and then apply that filter in a
                                % zero-phase compliant way.
                                participantAverages{ThisBin}(k,ThisChan,:) = filtfilt(dFilter,temp(ThisChan,:));
                                % figure out how many non-nan epochs there
                                % were.
                                validEpochs = length(squeeze(mean(GoodTrials(ThisBin).data(keyChanIdx,:,:), [2,1],'omitnan')));
                                epochCounts{ThisBin}(k) = validEpochs;
                            else
                                participantAverages{ThisBin}(k,ThisChan,:) = filtfilt(dFilter,temp(ThisChan,:));
                                validEpochs = length(squeeze(mean(GoodTrials(ThisBin).data(keyChanIdx,:,:), [2,1],'omitnan')));
                                epochCounts{ThisBin}(k) = validEpochs;
                            end
                            
                        else % no additional filters. Just take the data as they are.
                            participantAverages{ThisBin}(k,ThisChan,:) = temp(ThisChan,:);
                            %
                            validEpochs = length(squeeze(mean(GoodTrials(ThisBin).data(keyChanIdx,:,:), [2,1],'omitnan')));
                            epochCounts{ThisBin}(k) = validEpochs;
                        end
                        
                        % apply a time correction, if requested.
                        if ~isempty(sampleCorrection)
                            % time correction sought...
                            if sampleCorrection(ThisBin) == 0
                                % but not during this bin. skip.
                            else
                                % yep, need to adjust this paricular bin.
                                shft = abs(sampleCorrection(ThisBin));
                                tmp = squeeze(participantAverages{ThisBin}(k,ThisChan,:))';
                                if sampleCorrection(ThisBin) < 0
                                    % bring this bin forward in time.
                                    tmp = [tmp(shft+1:end),zeros(1,shft)];
                                else % push this bin back in time.
                                    tmp = [zeros(1,shft), tmp(1:end-shft)];
                                end
                                participantAverages{ThisBin}(k,ThisChan,:) = tmp;
                            end % of bin by bin check
                        end % sampleCorrection check.
                    end % of channel by channel loop.
                    
                    display(['Processing SUB ' SUB{k} ' Bin ' num2str(ThisBin)]);
                else
                    display(['Skipping SUB ' SUB{k} ' Bin ' num2str(ThisBin) '. Too few epochs.']);
                end
            end
        end % of skipping empty data sets
    end % of bin by bin loop.
end % of subject by subject loop


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
if ~exist('ERP_GrandAverages', 'dir')
    mkdir('ERP_GrandAverages');
end

keyPeriod = (times > wholeEpoch(1)/1000 & times < wholeEpoch(2)/1000);

% participantAverages.
% structure: a n-element cell array (cells are bins). Each bin contains a
% 3d structure: PIDs by chans by samples.
% average across channels (if more than one key channel chosen).
for ThisBin = 1:NoOfBins
    % pool the channels (if you have a multi channel montage).
    temp = squeeze(nanmean(participantAverages{ThisBin}(:,keyChanIdx,:),2));
    % limit to the relevant epoch time period.
    tempForEval{ThisBin} = temp; % will use for evaluation of key window two steps down.
    tempForOutput{ThisBin} = temp(:,keyPeriod);
end

% save that raw data (i.e. tempForOutput)
rawdataFilename = [pwd filesep 'ERP_GrandAverages' filesep 'RawOutput_PIDbySamples.xlsx'];
for ThisBin = 1:NoOfBins
    % declare a useful tab name
    if isempty(binContrast)
        tabname = ['Bin' num2str(ThisBin)];
    else
        tabname = 'DifferenceWave';
    end
    % write the data
    writematrix(tempForOutput{ThisBin}, rawdataFilename, 'Sheet', tabname, 'Range','B3');
    % write the row headers
    writecell(DataConfig.SUB', rawdataFilename, 'Sheet', tabname, 'Range', 'A3');
    writecell({'Times'}, rawdataFilename, 'Sheet', tabname, 'Range', 'A1');
    writecell({'PID'}, rawdataFilename, 'Sheet', tabname, 'Range', 'A2');
    % write the column headers
    writematrix(times(keyPeriod), rawdataFilename, 'Sheet', tabname, 'Range','B1');
end


% when do we want to quantify as the "key" period?
measurePeriod = (times > measureWindow (1)/1000 & times < measureWindow(2)/1000);

% and now do the same, but average across measurement window and outputting
% a CSV with PID and mean value in window.

% tempForEval is the variable we need to use now.
% structure: a n-element cell array (cells are bins). Each bin contains a
% 2d structure: PIDs by samples (measured at the declared chan).

processedFilename = [pwd filesep 'ERP_GrandAverages' filesep 'Processed_meanVoltagePerPerson.xlsx'];
for ThisBin = 1:NoOfBins
    % go bin by bin and take average of all samples within measurePeriod.
    tempForEval{ThisBin} = nanmean(tempForEval{ThisBin}(:,measurePeriod),2);
    % declare a useful tab name
    if isempty(binContrast)
        tabname = ['Bin' num2str(ThisBin)];
    else
        tabname = 'DifferenceWave';
    end
    % write that as output.
    writematrix(tempForEval{ThisBin}, processedFilename, 'Sheet', tabname, 'Range','B3');
    % write the epoch counts
    writematrix(epochCounts{ThisBin}', processedFilename, 'Sheet', tabname, 'Range','C3');
    % write the row headers
    writecell(DataConfig.SUB', processedFilename, 'Sheet', tabname, 'Range', 'A3');
    writecell({'PID'}, processedFilename, 'Sheet', tabname, 'Range', 'A2');
    % write the column headers
    writecell({'MeanVoltage'}, processedFilename, 'Sheet', tabname, 'Range', 'B2');
    writecell({'NoOfEpochs'}, processedFilename, 'Sheet', tabname, 'Range', 'C2');
end % of bin by bin loop.


%% and now start drawing.
% participantAverages.
% structure: participants, chans, samples, by bins.

% we want global mean and SEM over the keyPeriod time window.
% and we want separate figures for each bin.

for ThisBin = 1:NoOfBins
    timesToPlot = times(keyPeriod);
    % report grand average(across PIDs) for nominated chans in this
    % bin.
    % data in PID by samples format, so average across PIDs.
    meansToPlot = nanmean(tempForOutput{ThisBin},1);
    SEMs = std(tempForOutput{ThisBin},0,1, 'omitnan');
    % now should be: PIDs by samples
    minToPlot = meansToPlot - SEMs;
    maxToPlot = meansToPlot + SEMs;
    
    % find peak within measurement window?
    % PIDs by chans by samples.
    findingPeaks = squeeze(nanmean(participantAverages{ThisBin} (:,keyChanIdx,measurePeriod), [2, 1]));
    [posPeak, posPeakIdx] = max(findingPeaks);
    [negPeak, negPeakIdx] = min(findingPeaks);
    timesToMeasure = times(measurePeriod);
    posPeakTime = timesToMeasure(posPeakIdx);
    negPeakTime = timesToMeasure(negPeakIdx);
    posPeakText = ['Max Val = ' num2str(round(posPeak,3)) 'uV at ' ...
        num2str(round(posPeakTime,2)) 's' ];
    negPeakText = ['Min Val = ' num2str(round(negPeak,3)) 'uV at ' ...
        num2str(round(negPeakTime,2)) 's' ];
    
    % start drawing
    figure;
    hold on
    xline(0, '-k'); % show time zero
    xline(measureWindow(1)/1000, ':k'); % show start of eval period.
    xline(measureWindow(2)/1000, ':k'); % show end of eval period.
    if showIndividualTraces == 1
        for ThisPID = 1:size(tempForOutput{ThisBin},1)
            line(timesToPlot, tempForOutput{ThisBin}(ThisPID,:), 'LineStyle', ':', 'Color', 'k', 'LineWidth', 0.5);
        end
    end
    line(timesToPlot, meansToPlot, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 2);
    line(timesToPlot, minToPlot, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 2);
    line(timesToPlot, maxToPlot, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 2);
    hold off
    if isempty(binContrast)
        title(['Global Mean And SEMs Of Bin: ' num2str(ThisBin) ' at channel(s): ' string(strjoin(keyChans))]);
    else
        title(['Global Mean And SEMs Of Difference Wave at channel(s): ' string(strjoin(keyChans))]);
    end
    if showPeakInfo == 1
        % provide some info overlays about peaks.
        text(posPeakTime, posPeak, ['\leftarrow ' posPeakText], 'Color','green','FontSize',10);
        text(negPeakTime, negPeak, ['\leftarrow ' negPeakText], 'Color','green','FontSize',10);
    end
    ylabel('Voltage(microvolts)');
    xlabel('Time(s)');
    y_cap = 2*max(abs(meansToPlot));
    if isempty(fixedYlim)
      ylim([-1*y_cap, y_cap]);  
    else
      ylim(fixedYlim);  
    end
    f = gcf;
    f.Units = 'inches';
    f.OuterPosition = [0.5 0.5 5.5 5.5]; % make the figure 5 inches in size.
    fig_filename = ['ERP_GrandAverages' filesep 'Bin_' num2str(ThisBin) '_GrandAverage.png'];
    disp(['Saving ERP image ' fig_filename]);
    exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
    close(gcf);
    
    % and now do a topoplot per bin.
    
    % ok, need to pull mean voltage over measurement window, per
    % channel. And that needs to be averaged across participants too.
    
    % participantAverages
    % structure: a n-element cell array (cells are bins). Each bin contains a
    % 3d structure: PIDs by chans by samples.
    % sequentially averages by saampels, then across PIDs. But only in
    % measurement window.
    dataToPlot = nanmean(participantAverages{ThisBin} (:,:,measurePeriod), [3, 1]);
    % limit the plotting to scalp channels.
    dataToPlot = dataToPlot(1:DataConfig.TotalChannels{1});
    % now should be: PIDs by chans
    
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
            topoplot(dataToPlot, chanlocs(1:DataConfig.TotalChannels{1}), 'electrodes' ,'ptslabels');
        else
            % you just want a standard image. No info overlays.
            topoplot(dataToPlot, chanlocs(1:DataConfig.TotalChannels{1}), 'style', 'fill', 'numcontour', 3);

            % topoplot(dataToPlot, chanlocs);
        end
        colorbar;
        colormap(colScheme);
        if isempty(binContrast)
            title(['Topoplot of Bin: ' num2str(ThisBin) 'during measurement window']);
        else
            title(['Topoplot of difference wave during measurement window']);
        end
        colourCap = 1.5*max(abs(dataToPlot));
        
        if isempty(fixedYlim)
            caxis([-1*colourCap, colourCap]);
        else
            caxis(fixedYlim);
        end
        
        if showPeakInfo == 1
            % add text related to min/max
            annotation('textbox', [0, 0.15, 0, 0], 'string', minText, 'FitBoxToText','on');
            annotation('textbox', [0, 0.05, 0, 0], 'string', maxText, 'FitBoxToText','on');
        end
        f = gcf;
        f.Units = 'inches';
        f.OuterPosition = [0.5 0.5 5.5 5.5]; % make the figure 5 inches in size.
        fig_filename = ['ERP_GrandAverages' filesep 'Bin_' num2str(ThisBin) '_GrandTopoplot.png'];
        disp(['Saving topoimage image ' fig_filename]);
        exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
        close(gcf);
    end
end % of cycling through bins.

% draw one more figure that captures all bins using default matlab
% settings. Good for comparing across bins, but only if there are multiple
% bins.


% a set of six colours (which then repeat, if necessary)
lineColours = {[1, 0, 0], [0, 1, 0], [0, 0, 1], ...
    [1, 1, 0], [0, 1, 1], [1, 0, 1]};
shadeColours = {[0.8, 0.2, 0.2], [0.2, 0.8, 0.2], [0.2, 0.2, 0.8], ...
    [0.8, 0.8, 0.2], [0.2, 0.8, 0.8], [0.8, 0.2, 0.8]};

if isempty(binContrast) & NoOfBins > 1
    dataToPlot = zeros(NoOfBins,length(timesToPlot));
    figure;
    hold on
    for ThisBin = 1:NoOfBins
        mean_vals = mean(tempForOutput{ThisBin},1, 'omitnan');
        BinLbls{ThisBin} = ['Bin ' num2str(ThisBin)];
        % find PID count.
        validPIDs = sum(~isnan(tempForOutput{ThisBin}),2);
        PIDcount = length(validPIDs(validPIDs > 0));
        %
        plus_sem = mean_vals + ...
            std(tempForOutput{ThisBin},[],1, 'omitnan')./sqrt(PIDcount-1);
        minus_sem = mean_vals - ...
            std(tempForOutput{ThisBin},[],1,'omitnan')./sqrt(PIDcount-1);
        % calculate the shaded region
        x2 = [timesToPlot, fliplr(timesToPlot)];
        inBetween = [plus_sem, fliplr(minus_sem)];
        
        % make an adjustment to accomodate > 6 conditions.
        thisLine = mod(ThisBin,6);
        if thisLine == 0
            thisLine = 6;
        end
        
        fill(x2, inBetween,shadeColours{thisLine}, 'FaceAlpha',0.5, 'LineStyle', 'none');
        
        % and now draw the mean line plot
        line(timesToPlot, mean_vals, ...
            'Color', lineColours{thisLine}, 'LineWidth', 3, 'LineStyle', '-' );
    end
    
    % and sort out some other peripherals
    % title('AllBins');
    % legend(BinLbls);
    xline(0, 'LineWidth', 1.5,'LineStyle', ':'); % show time zero
    ylabel('Voltage(microvolts)', 'Fontsize', 16);
    xlabel('Time(s)', 'Fontsize', 16);
    ax = gca;
    ax.FontSize = 16;
    y_cap = 2*max(abs(dataToPlot(:)));
    if isempty(fixedYlim)
        ylim([-1*y_cap, y_cap]);
    else
        ylim(fixedYlim);
    end
    %
    f = gcf;
    f.Units = 'inches';
    f.OuterPosition = [0.5 0.5 5.5 5.5]; % make the figure 5 inches in size.
    fig_filename = ['ERP_GrandAverages' filesep 'AllBinsTogether.png'];
    disp('Saving all bin waveform');
    exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
    close(gcf);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Danielle specific drawing script.
% a set of six colours (which then repeat, if necessary)

% which bins to contrast with each other?
outFigs = {[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [2, 3, 8, 9]};
condlbls = {'AudStand', 'AudDevRew', 'AudDevCtrl', ...
    'VisStand', 'VisDevRew', 'VisDevCtrl', ...
    'BaseAudStand', 'BaseAudDevRew', 'BaseAudDevCtrl', ...
    'BaseVisStand', 'BaseVisDevRew', 'BaseVisDevCtrl'};

for thisFig = 1:length(outFigs)
    
    lineColours = {[1, 0, 0], [0, 1, 0], [0, 0, 1], ...
        [1, 1, 0], [0, 1, 1], [1, 0, 1]};
    shadeColours = {[0.8, 0.2, 0.2], [0.2, 0.8, 0.2], [0.2, 0.2, 0.8], ...
        [0.8, 0.8, 0.2], [0.2, 0.8, 0.8], [0.8, 0.2, 0.8]};
    
    if isempty(binContrast) & NoOfBins > 1
        dataToPlot = zeros(NoOfBins,length(timesToPlot));
        figure;
        hold on
        binsToPlot = outFigs{thisFig};
        clear BinLbls;
        
        % shade in the measurement region.
        xShade = [measureWindow(1)/1000, measureWindow(1)/1000, ...
            measureWindow(2)/1000, measureWindow(2)/1000];
        yShade = [fixedYlim(1), fixedYlim(2), fixedYlim(2), fixedYlim(1)];
        fill(xShade, yShade, [0.5, 0.5, 0.5], 'FaceAlpha',0.5, 'LineStyle', 'none');
        BinLbls{1} = 'MeasureWindow';
            
        for ThisBin = 1:length(binsToPlot)
            mean_vals = mean(tempForOutput{binsToPlot(ThisBin)},1, 'omitnan');
            BinLbls{ThisBin*2} = condlbls{binsToPlot(ThisBin)};
            BinLbls{ThisBin*2+1} = condlbls{binsToPlot(ThisBin)};
            % find PID count.
            validPIDs = sum(~isnan(tempForOutput{binsToPlot(ThisBin)}),2);
            PIDcount = length(validPIDs(validPIDs > 0));
            %
            plus_sem = mean_vals + ...
                std(tempForOutput{binsToPlot(ThisBin)},[],1, 'omitnan')./sqrt(PIDcount-1);
            minus_sem = mean_vals - ...
                std(tempForOutput{binsToPlot(ThisBin)},[],1,'omitnan')./sqrt(PIDcount-1);
            % calculate the shaded region
            x2 = [timesToPlot, fliplr(timesToPlot)];
            inBetween = [plus_sem, fliplr(minus_sem)];
            
            % make an adjustment to accomodate > 6 conditions.
            thisLine = mod(ThisBin,6);
            if thisLine == 0
                thisLine = 6;
            end
            
            fill(x2, inBetween,shadeColours{thisLine}, 'FaceAlpha',0.5, 'LineStyle', 'none');
            
            % and now draw the mean line plot
            line(timesToPlot, mean_vals, ...
                'Color', lineColours{thisLine}, 'LineWidth', 3, 'LineStyle', '-' );
        end
        
        % and sort out some other peripherals
        % title('AllBins');
        lgd = legend(BinLbls, 'Location', 'northoutside');
        lgd.FontSize = 8;
        xline(0, 'LineWidth', 1.5,'LineStyle', ':'); % show time zero
        ylabel('Voltage(microvolts)', 'Fontsize', 16);
        xlabel('Time(s)', 'Fontsize', 16);
        ax = gca;
        ax.FontSize = 16;
        y_cap = 2*max(abs(dataToPlot(:)));
        if isempty(fixedYlim)
            ylim([-1*y_cap, y_cap]);
        else
            ylim(fixedYlim);
        end
        %
        f = gcf;
        f.Units = 'inches';
        f.OuterPosition = [0.5 0.5 7.5 7.5]; % make the figure 7 inches in size.
        fig_filename = ['ERP_GrandAverages' filesep 'Bins_' ...
            num2str(outFigs{thisFig}) '.png'];  
        disp('Saving all bin waveform');
        exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
        close(gcf);
    end
    
end

