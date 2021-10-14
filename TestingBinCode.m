%% set up a config file with relevant SUBS information
% let's loop through all relevant files and draw the ERG1 figures.
Current_File_Path = pwd;
addpath('Functions');
ConfigFileName = 'WIMR_Config_testing';
ConfigFilePath = [Current_File_Path filesep 'SupportingDocs' filesep ConfigFileName '.xlsx'];
Options = detectImportOptions(ConfigFilePath);

for k = 1:numel(Options.VariableTypes)
    Options.VariableTypes{k} = 'char';
end
DataConfig = table2struct(readtable(ConfigFilePath, Options));
DataConfig = adjustConfigData(DataConfig);

% and open eeglab to access those functions
eeglab;

%% loop through SUBS and draw relelvant figures.
SUB = DataConfig.SUB;
for k = 1:length(SUB)
    
    % just working space to check how data is binned in X7 code.
    testFolder = [fileparts(pwd) filesep SUB{k}];
    testFile_EEG = [SUB{k} '_ds_PREP_ica_corr_cbip_elist_bins_epoch_ar.set'];
    testFile_mat = [SUB{k} '_ARcorrectedBins.mat'];

    % channels to run test on.
    FCz = 47;
    Cz = 48;
    Erg1 = 72;
    Erg2 = 73;
    Bin = 1;
    avFCz = [ 11, 38, 46, 47, 48];
    Mast1= 65;
    Mast2 = 66;
    
    % load the outputted data in matlab format.
    load([testFolder filesep testFile_mat]);
    % creates structure "GoodTrials"
    % load the raw data in EEGlab format.
    EEG = pop_loadset(testFile_EEG, testFolder);
    
    % generate an x-axis for plotting.
    times = EEG.times;
    % times = ([1:length(dataToPlot1)].*1/EEG.srate) - 0.2;
    
    keyPeriod = (times > -200 & times < 800);
    
    for Bin = 1:6
        if Bin > numel(GoodTrials) || isempty(GoodTrials(Bin).data) || isempty(EEG.data)
            % then nothing to analyze or draw. 
        else
            temp = squeeze(mean(GoodTrials(Bin).data(Mast1,keyPeriod,:),1));
            dataToPlot2 = mean(temp,2);
            temp = squeeze(mean(GoodTrials(Bin).data(Mast2,keyPeriod,:),1));
            dataToPlot3 = mean(temp,2);
            timesToPlot = times(keyPeriod);
            figure;
            hold on
            line(timesToPlot, dataToPlot2, 'Color', 'blue');
            line(timesToPlot, dataToPlot3, 'Color', 'black');
            hold off
            title(['Bin_' num2str(Bin)]);
            outFigureName = [testFolder filesep 'Figures' filesep ...
                'PID_' SUB{k} '_Bin_' num2str(Bin) '_compareMasts.png'];
            saveas(gcf,outFigureName);
            close(gcf);
            % and then do the same, but for Erg2
            temp = squeeze(mean(GoodTrials(Bin).data([Mast1, Mast2],keyPeriod,:),1));
            dataToPlot2 = mean(temp,2);
            temp = squeeze(mean(EEG.data([Mast1, Mast2],keyPeriod,:),1));
            dataToPlot3 = mean(temp,2);
            timesToPlot = times(keyPeriod);
            figure;
            hold on
            line(timesToPlot, dataToPlot2, 'Color', 'blue');
            line(timesToPlot, dataToPlot3, 'Color', 'black');
            hold off
            title(['Bin_' num2str(Bin)]);
            outFigureName = [testFolder filesep 'Figures' filesep ...
                'PID_' SUB{k} '_Bin_' num2str(Bin) '_compareAverageMasts.png'];
            saveas(gcf,outFigureName);
            close(gcf);
        end % of skipping out for empty data sets.
    end
    
end % of subject by subject loop


