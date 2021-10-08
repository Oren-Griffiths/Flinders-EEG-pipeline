% simple little function that pulls out the cleaned data per bin in the
% format EEG = (channels, samples, epochs). there is one dataset per bin.
% T

function GoodData = X7_AverageERPs

% initialize
global DataConfig
close all; 
clearvars -except DataConfig;
SUB = DataConfig.SUB;

%Location of the main study directory
DIR = fileparts(fileparts(mfilename('fullpath'))); 

%Location of the folder that contains this script and any associated processing files
Current_File_Path = fileparts(mfilename('fullpath'));


%Open EEGLAB and ERPLAB Toolboxes  
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Loop through each subject listed in SUB
for i = 1:length(SUB)
tic;

DIR = fileparts(fileparts(mfilename('fullpath')));
Subject_Path = [DIR filesep SUB{i} filesep];

EEG = pop_loadset( 'filename', [ SUB{i} ...
    '_ds_reref_ucbip_hpfilt_interp_event_ica_corr_cbip_elist_bins_epoch_ar.set'], 'filepath', Subject_Path);

% NEXT pull out different output matrices per bin...
% Do this by looking in EEG.epoch(ThisEpoch).eventbini
% The entries will either be -1 (not a bin label) or a positive integer
% referring to bin that event belongs to: 1,2,3, etc. Stored as double.

% so search for first non "-1" value. If reach end, not in a bin. If in a
% bin, you'll hit a positive integer pretty quickly.

% Then you know which bucket to put that trial in. Probably have to create
% a numeric index to add to EEG.data (channels, samples, epochs) to make it
% (channels, sample, epochs) plus a (bins) vector of same lengths epochs
% and then you can just cut the bits out that you need after "AR-failed"
% trials are removed.

% Finds all the non AR-rejected trials.
% start by initializing some output vectors.
FlaggedEpochs = ones(numel(EEG.epoch),1);
EpochBins = zeros(numel(EEG.epoch),1);

for ThisTrial = 1:numel(EEG.epoch)
    % loop through every epoch
    FlagsThisTrial_cell = EEG.epoch(ThisTrial).eventflag; % collect flags in cell
    BinThisCell_cell = EEG.epoch(ThisTrial).eventbini;
    % empty flags (0) are double, but AR-flagged items are unit16
    % so can't vectorize, must loop. So let's loop.
    
    % Loop through the possible flag values and bin values.
    
    for ThisCell = 1:numel(FlagsThisTrial_cell)
        if double(any(FlagsThisTrial_cell{ThisCell})) > 0
            FlaggedEpochs(ThisTrial) = 0;
        end
        if EpochBins(ThisTrial) > 0 % Bin number already found, so skip this.
        else
            if BinThisCell_cell{ThisCell} > 0
                EpochBins(ThisTrial) = double(BinThisCell_cell{ThisCell});
            end
        end
    end
    
    
end

% put all the non-rejected trials in one matrix and reports on time taken.
TotalGoodTrials = EEG.data(:,:,logical(FlaggedEpochs));
GoodEpochBins = EpochBins(logical(FlaggedEpochs));
GoodTrials = {};
BinIDs = unique(GoodEpochBins)
for ThisBin = BinIDs'
    disp(['Parsing bin number ' num2str(ThisBin)]);
    GoodTrials{ThisBin} = EEG.data(:,:,GoodEpochBins == ThisBin);
end
PropnFlagged = nnz(~FlaggedEpochs)/length(FlaggedEpochs);
TimeTaken = toc;

disp(['Number of flagged trials is ' num2str(nnz(~FlaggedEpochs))]);
disp(['Propn of flagged trials is ' num2str(PropnFlagged)]);
disp(['Time taken to do this is ' num2str(toc) 'seconds']);
end % subject by subject loop
end % end function.

