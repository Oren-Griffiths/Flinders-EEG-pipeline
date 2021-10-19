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

% what time period (ms) do you want visualized. This will obviously break
% if you choose an area large than the epoch declared in DataConfig, so
% choose sensibly.
wholeEpoch = [-200, 800];

% choose a time period you want to take the average across. Measured in ms.
measureWindow = [50, 150];

% if a person has less than X clean epochs are AR-rejection, then remove
% them from the averaging process. Set X. Put empty '[]' to ignore min.
minEpochs = 30;

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
end

chanlocs = GoodTrials(1).chanlocs; % we need this structure too.

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
% structure: participants, chans, samples, by bins.
participantAverages= NaN(length(DataConfig.SUB),NoOfChans,LengthOfEpoch, NoOfBins);

% generate an x-axis for plotting. Measurement is now in seconds.
times = ([1:LengthOfEpoch].*1/srate) + (DataConfig.EpochMin{1}/1000);

%% loop through SUBS and gather *per participant* averages by averaging across epochs (within bins).
SUB = DataConfig.SUB;
for k = 1:length(SUB)
    
    % just working space to check how data is binned in X7 code.
    testFolder= [fileparts(pwd) filesep SUB{k}];
    testFile_mat = [SUB{k} '_ARcorrectedBins.mat'];
    
    % load the outputted data in matlab format.
    load([testFolder filesep testFile_mat]);
    
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
                        participantAverages(k,ThisChan,:,ThisBin) =  temp(ThisChan,:);
                    end
                    display(['Processing SUB ' SUB{k} ' Bin ' num2str(ThisBin)]);  
                else
                    display(['Skipping SUB ' SUB{k} ' Bin ' num2str(ThisBin) '. Too few epochs.']);  
                end
                
            end
        end % of skipping empty data sets
    end % of bin by bin loop.
end % of subject by subject loop

%% and now start to calculate the output the data needed.
if ~exist('ERP_GrandAverages', 'dir')
    mkdir('ERP_GrandAverages');
end

keyPeriod = (times > wholeEpoch(1)/1000 & times < wholeEpoch(2)/1000);
% get channel numbers
% that is, convert channel indices from channel names.
% have done this assuming Bin 1 is the same as for all others, but possible
% to adapt to do this independently per bin if montage changes across
% bins(?)
for ThisChan = 1:length(keyChans)
    keyChanIdx(ThisChan) = find(strcmp({GoodTrials(1).chanlocs.labels}, keyChans{ThisChan})==1);
end

% participantAverages.
% structure: participants, chans, samples, by bins.
% average across channels (if more than one key channel chosen).
temp = squeeze(nanmean(participantAverages(:,keyChanIdx,:,:),2));
if length(DataConfig.SUB) == 1
    % if only one participant, that dimension will be "squeezed" too.
    %  turn this into 3D structure for consistency.
    temp2 = ones(1,size(temp,1), size(temp, 2));
    temp2(1,:,:) = temp(:,:);
    temp = temp2;
end

% now structure is: participants, samples, by bins.
% limit to the relevant epoch.
tempForOutput = temp(:,keyPeriod,:);

% save that raw data
rawOutput = [pwd filesep 'ERP_GrandAverages' filesep 'ChosenChans_PIDbySamplesByBins.mat'];
save(rawOutput, 'tempForOutput');

% when do we want to quantify as the "key" period?
measurePeriod = (times > measureWindow (1)/1000 & times < measureWindow(2)/1000);

% and now do the same, but average across measurement window and outputting
% a CSV with PID and mean value in window.
PIDxSampxBin = temp;
for ThisBin = 1:NoOfBins
    for k = 1:length(SUB)
        meanByBin{k+1,1} = SUB{k};
        meanByBin{k+1,ThisBin+1} = nanmean(PIDxSampxBin(k, measurePeriod, ThisBin),2);
    end % of PID by PID loop.
    % add a header per bin
    meanByBin{1,ThisBin+1} = ['Bin' num2str(ThisBin)];
end % of bin by bin loop.
% and just one more header
meanByBin{1,1} = 'PID';

% and output that value
meansOutput = [pwd filesep 'ERP_GrandAverages' filesep 'ValuesInMeasureWindow.xlsx'];
writecell(meanByBin, meansOutput);

%% and now start drawing.
% participantAverages.
% structure: participants, chans, samples, by bins.

% we want global mean and SEM over the keyPeriod time window.
% and we want separate figures for each bin.
inputSize = size(participantAverages);
temp = [];
if min(inputSize) > 1
    for ThisBin = 1:NoOfBins
        timesToPlot = times(keyPeriod);
        % average across nominated chans (if more than one).
        temp = squeeze(nanmean(participantAverages(:,keyChanIdx,keyPeriod,ThisBin),2));
        % now should be: PIDs by samples
        meansToPlot = nanmean(temp,1);
        SEMs = std(temp,0,1, 'omitnan');
        minToPlot = meansToPlot - SEMs;
        maxToPlot = meansToPlot + SEMs;
        figure;
        hold on
        xline(0, ':k'); % show time zero
        xline(measureWindow(1)/1000, ':k'); % show start of eval period.
        xline(measureWindow(2)/1000, ':k'); % show end of eval period.
        for ThisPID = 1:size(temp,1)
            line(timesToPlot, temp(ThisPID,:), 'LineStyle', '-', 'Color', 'k', 'LineWidth', 0.5);
        end
        
        line(timesToPlot, meansToPlot, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 3);
        line(timesToPlot, minToPlot, 'LineStyle', ':', 'Color', 'r', 'LineWidth', 2);
        line(timesToPlot, maxToPlot, 'LineStyle', ':', 'Color', 'r', 'LineWidth', 2);
        hold off
        title(['Global Mean And SEMs Of Bin: ' num2str(ThisBin) ' at channel(s): ' string(strjoin(keyChans))]);
        ylabel('Voltage(microvolts)');
        xlabel('Time(s)');
        % change some values and save.
        f = gcf;
        f.Units = 'inches';
        f.OuterPosition = [0.5 0.5 5.5 5.5]; % make the figure 5 inches in size.
        fig_filename = ['ERP_GrandAverages' filesep 'Bin_' num2str(ThisBin) '_GrandAverage.png'];
        disp(['Saving ERP image ' fig_filename]);
        exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
        close(gcf);
        
        % and now do a topoplot per bin.
        measurePeriod;
        temp = [];
        temp = squeeze(nanmean(participantAverages(:,1:DataConfig.TotalChannels{1},measureWindow,ThisBin),3));
        dataToPlot = squeeze(nanmean(temp,1));
        % now should be: PIDs by chans
        
        if sum(isnan(dataToPlot)) == length(dataToPlot)
            % don't actually draw the image if there aren't any data.
        else 
            figure;
            topoplot(dataToPlot, chanlocs);
            colorbar;
            title(['Topoplot of Bin: ' num2str(ThisBin) 'during measurement window']);
            f = gcf;
            f.Units = 'inches';
            f.OuterPosition = [0.5 0.5 5.5 5.5]; % make the figure 5 inches in size.
            fig_filename = ['ERP_GrandAverages' filesep 'Bin_' num2str(ThisBin) '_GrandTopoplot.png'];
            disp(['Saving topoimage image ' fig_filename]);
            exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
            close(gcf);
        end
    end
else
    disp('Either participants, chans, samples, or bins is size 1, so cannot draw figures.');
end

%% ok, last stage is to do the relevant topoplots


