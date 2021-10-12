% just working space to check how data is binned in X7 code. 
testFolder= [fileparts(pwd) filesep '810'];
testFile_EEG = ['810_ds_PREP_ica_corr_cbip_elist_bins_epoch_ar.set'];
testFile_mat = ['810_ARcorrectedBins.mat'];

% channels to run test on.
FCz = 47; 
Cz = 48;
Erg1 = 72; 
Erg2 = 73; 
Bin = 1;

% load the outputted data in matlab format. 
load([testFolder filesep testFile_mat]);
% creates structure "GoodTrials"
% load the raw data in EEGlab format.
EEG= pop_loadset(testFile_EEG, testFolder);

% plot bin X at FCz for both files.
dataToPlot1 = squeeze(mean(GoodTrials(Bin).data(FCz,:,:),3))';
% generate an x-axis for plotting.
times = ([1:length(dataToPlot1)].*1/256) - 0.2;
figure; plot(times,dataToPlot1);
title('Extracted_Bin1data');

RawGoodTrials = GoodTrials;


for Bin = 1:6
dataToPlot1 = squeeze(mean(RawGoodTrials(Bin).data(FCz,:,:),3))';
dataToPlot2 = squeeze(mean(GoodTrials(Bin).data(FCz,:,:),3))';
dataToPlot3 = squeeze(mean(EEG.data(FCz,:,:),3))';
figure;
hold on
line(times, dataToPlot1, 'Color', 'r');
line(times, dataToPlot2, 'Color', 'b');
line(times, dataToPlot3, 'Color', 'black'); 
hold off
end
%
% dataToPlot2 = squeeze(mean(EEG.data(FCz,:,:),3));
% times = ([1:length(dataToPlot2)].*1/256) - 0.2;
% figure; plot(times,dataToPlot2);
% title('EEGlab_Bin1data');

