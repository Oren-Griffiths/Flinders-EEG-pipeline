% simple little function that pulls out the cleaned data per bin in the
% format EEG = (channels, samples, epochs). there is one dataset per bin.

% updated: 13.09.21
% primary update is to make the function output an image per bin. That
% image shows the data (that passes AR filter) plotted epoch-by-epoch.
% Allows user to see if a few errant epochs are driving the mean. Currently
% plots channel Cz or, if specified in ConfigFile, will plot the channel
% chose in DataConfig.KeyChannel{1}.

function X7_ExtractEpochedData_p(DataConfig,SUB, imageType)

% initialize
close all;
clearvars -except DataConfig SUB imageType;

% Location of the main study directory
DIR = fileparts(pwd);

% location of preprocessing files.
Current_File_Path = pwd;


%Open EEGLAB and ERPLAB Toolboxes
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Loop through each subject listed in SUB
for i = 1:length(SUB)
    tic;
    
    Subject_Path = [DIR filesep SUB{i} filesep];
    
    EEG = pop_loadset( 'filename', [ SUB{i} ...
        '_ds_PREP_ica_corr_cbip_elist_bins_epoch_ar.set'], 'filepath', Subject_Path);
    
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
        % can't vectorize, must loop. So let's loop.
        
        % Loop through the possible flag values and bin values.
        for ThisCell = 1:numel(FlagsThisTrial_cell)
            if double(any(FlagsThisTrial_cell{ThisCell})) > 0
                % exclude this trial as there's an artefact.
                FlaggedEpochs(ThisTrial) = 0;
            end
            if EpochBins(ThisTrial) > 0
                % Bin number already found, so skip this check.
            else % look for a bin number.
                if BinThisCell_cell{ThisCell} > 0
                    EpochBins(ThisTrial) = double(BinThisCell_cell{ThisCell});
                end
            end
        end
        
        
    end % of loop through every epoch (to remove flagged epochs)
    
    % put all the non-rejected trials in one matrix and reports on time taken.
    TotalGoodTrials = EEG.data(:,:,logical(FlaggedEpochs));
    GoodEpochBins = EpochBins(logical(FlaggedEpochs));
    BinIDs = unique(GoodEpochBins);
    GoodTrials = [];
    if length(BinIDs) > 0
        for idx = 1:length(BinIDs')
            ThisBin = BinIDs(idx);
            disp(['Parsing bin number ' num2str(ThisBin)]);
            % find artefact-free bins that I want to scoop up.
            epochsToTake = GoodEpochBins == ThisBin;
            if sum(EpochBins == ThisBin) > 0
                rejRate = 1 - (sum(GoodEpochBins == ThisBin)/ ...
                    sum(EpochBins == ThisBin));
            else
                rejRate = 0;
            end
            % then scoop them up.
            GoodTrials(ThisBin).data = double(TotalGoodTrials(:,:,epochsToTake));
            GoodTrials(ThisBin).rejRate = rejRate;
            GoodTrials(ThisBin).ID = ThisBin;
            GoodTrials(ThisBin).chanlocs = EEG.chanlocs;
            GoodTrials(ThisBin).srate = EEG.srate;
        end
    else % no passing trials, so just note which events occured with 0 data
        for idx = 1:length(EpochBins)
            ThisBin = EpochBins(idx);
            GoodTrials(ThisBin).data = [];
            GoodTrials(ThisBin).rejRate = 1;
            GoodTrials(ThisBin).ID = ThisBin;
            GoodTrials(ThisBin).chanlocs = EEG.chanlocs;
            GoodTrials(ThisBin).srate = EEG.srate;
        end
    end
    PropnFlagged = nnz(~FlaggedEpochs)/length(FlaggedEpochs);
    TimeTaken = toc;
    
    disp(['Number of flagged trials is ' num2str(nnz(~FlaggedEpochs))]);
    disp(['Propn of flagged trials is ' num2str(PropnFlagged)]);
    disp(['Time taken to do this is ' num2str(toc) 'seconds']);
    
    saveDestination = [Subject_Path  SUB{i} '_ARcorrectedBins.mat'];
    save(saveDestination, 'GoodTrials');
    
    if strcmp(imageType, 'none')
        % don't draw any images. Just save output data in .mat format. 
    else
        
        % draw one figure with a grand average of Fz, Cz, Pz of the good bins.
        figure;
        
        NoOfBins = numel(GoodTrials);
        if NoOfBins == 0
            % don't draw an average figure.
        else
            labels = {'Fz', 'Cz', 'Pz'};
            means = [];
            for ThisBin = 1:NoOfBins
                if ~isempty(GoodTrials(ThisBin).data)
                    % calculate grand mean of Fz, Cz, Pz. Draw a picture.
                    switch DataConfig.TotalChannels{1}
                        case 32
                            Chans = [31, 32, 13];
                            LineDetails = {'-r', '-b', '-g'};
                        case 64
                            Chans = [38, 48, 31];
                            LineDetails = {'-r', '-b', '-g'};
                    end
                    
                    % start drawing.
                    for idx = 1:length(Chans)
                        ThisChan = Chans(idx);
                        times = EEG.times;
                        means(idx,:) = mean(GoodTrials(ThisBin).data(ThisChan,:,:),3);
                        subplot( round(sqrt(NoOfBins))+1, round(sqrt(NoOfBins))+1, ThisBin);
                    end
                    for idx = 1:length(Chans)
                        hold on
                        plot(times,means(idx,:),LineDetails{idx});
                        hold off
                    end
                    legend(labels);
                    title(num2str(ThisBin));
                    set(gca,'FontSize',18);
                    
                end % of skipped empty bins
            end % of bin by bin loop
            set(gcf,'PaperPositionMode','manual','PaperUnits','Inches','PaperSize',[25 25], 'PaperPosition', [0 0 25 25] );
            out_filename = [Subject_Path 'Figures\X7_',SUB{i} , '_BinGrandMeans.pdf'];
            saveas(gcf,out_filename,'pdf');
            close(gcf);
        end
        
        % draw a second set of figures for each trial at "KeyChannel" overlaid
        % on each other (to get a sense of noise in data set for that bins).
        NoOfBins = numel(GoodTrials);
        if NoOfBins > 0
            % need to make sure that one "key" channel is specified.
            if isfield(DataConfig, 'KeyChannel')
                % no need to declare KeyChannel as it's already
                % specified.
                KeyChannel = DataConfig.KeyChannel{1};
            else
                % older versions of "MasterFile" may not declare KeyChannel config.
                % so just use Cz for them.
                switch DataConfig.TotalChannels{1}
                    case 32
                        KeyChannel = 32;
                    case 64
                        KeyChannel = 48;
                end
            end
            
            % loop through bins and draw a figure per bin.
            for ThisBin = 1:NoOfBins
                if ~isempty(GoodTrials(ThisBin).data)
                    % loop through every epoch at assigned channel.
                    figure;
                    % find x axis
                    times = EEG.times;
                    % grab the right bin.
                    dataToPlot = GoodTrials(ThisBin).data;
                    % squeeze down to epochs by samples.
                    if size(size(dataToPlot))< 3
                        % only one epoch retained, so no "epochs" dimension.
                        dataToPlot = dataToPlot(KeyChannel,:);
                    else % else index key channel to get samples by epochs.
                        dataToPlot = squeeze(dataToPlot(KeyChannel,:,:))';
                    end
                    % start drawing.
                    hold on
                    for k = 1:size(dataToPlot,1)
                        plot(times,dataToPlot(k,:),'-b');
                    end
                    % and overlay the mean.
                    plot(times, mean(dataToPlot), 'Color', 'r', 'LineWidth', 3 );
                    hold off
                    title(num2str(ThisBin));
                    set(gca,'FontSize',18);
                    
                    % and save it as it's own file.
                    set(gcf,'PaperPositionMode','manual','PaperUnits','Inches','PaperSize',[25 25], 'PaperPosition', [0 0 25 25] );
                    out_filename = [Subject_Path 'Figures\X7_PID_', SUB{i} , '_Bin_' num2str(GoodTrials(ThisBin).ID) '_RawData.png'];
                    saveas(gcf,out_filename);
                    close(gcf);
                    
                end % of skipped empty bins
            end % of bin-by-bin loop.
        end % of image skipping if 'none' selected
    end % of skipped if no data at all.
    
    
end % subject by subject loop

end % end function.

