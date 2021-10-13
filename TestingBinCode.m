% just working space to check how data is binned in X7 code. 
testFolder= [fileparts(pwd) filesep '810'];
testFile_EEG = ['810_ds_PREP_ica_corr_cbip_elist_bins_epoch_ar.set'];
testFile_mat = ['810_ARcorrectedBins.mat'];

% open eeglab to access those functions
eeglab;

% channels to run test on.
FCz = 47; 
Cz = 48;
Erg1 = 72; 
Erg2 = 73; 
Bin = 1;
avFCz = [ 11, 38, 46, 47, 48];

% load the outputted data in matlab format. 
load([testFolder filesep testFile_mat]);
% creates structure "GoodTrials"
% load the raw data in EEGlab format.
EEG = pop_loadset(testFile_EEG, testFolder);

% generate an x-axis for plotting.
times = EEG.times; 
% times = ([1:length(dataToPlot1)].*1/EEG.srate) - 0.2;

keyPeriod = (times > -800 & times < 800);

for Bin = 1:6
    temp = squeeze(mean(GoodTrials(Bin).data(Erg1,keyPeriod,:),1));
    dataToPlot2 = mean(temp,2);
    temp = squeeze(mean(EEG.data(Erg1,keyPeriod,:),1));
    dataToPlot3 = mean(temp,2);
timesToPlot = times(keyPeriod);
figure;
hold on
line(timesToPlot, dataToPlot2, 'Color', 'blue');
line(timesToPlot, dataToPlot3, 'Color', 'black'); 
hold off
title(['Bin_' num2str(Bin)]);
saveas(gcf,['Bin_' num2str(Bin) '_Erg1.png']); 
end

