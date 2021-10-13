% Oren's preprocessing script, based on ERP CORE MMN (Luck, Kappenman) and
% core functions of the PrepPipeline toolbox.

% update 15.09.21 to correct some calls to DataConfig.PREP{1} and also to
% make sure ERG channels untouched by cleaning.

function Y1_preprocess_wPREP

% initialize
global DataConfig
close all;
clearvars -except DataConfig;
SUB = DataConfig.SUB;

    % find the directory one up from this file (for subject folders)
    DIR = fileparts(pwd);
    
    % find the directory that this file lives in.
    Current_File_Path = pwd;
    
    
    %% Key input and preprocess loop
    
    %Loop through each subject listed in SUB
    for i = 1:length(SUB)
        
        %Open EEGLAB and ERPLAB Toolboxes
        %Acts as initialization of the relevant variables per participant.
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        
        %Define subject path based on study directory and subject ID of current subject
        Subject_Path = [DIR filesep num2str(SUB{i}) filesep];
        
        % record where we're up to in case of crash.
        DataConfig.CurrentSUB = SUB{i};
        
        % this will almost always be first process called, but just in case
        % it's not let's build in the check to figure out what file to open.
        
        if isempty(DataConfig.LastProcess)
            % this is the first process for this file, so open the raw data.
            if exist([Subject_Path  SUB{i} '_raw.set']) == 2
                % you've already imported this data set so just start
                % there.
                EEG = pop_loadset( 'filename',[SUB{i} '_raw.set'], 'filepath', Subject_Path);
            else
                switch DataConfig.RawFileType{1}
                    case '.bdf'
                        FileToOpen = [Subject_Path  SUB{i} '.bdf'];
                        if isfile(FileToOpen)
                            EEG = pop_biosig(FileToOpen, 'ref', DataConfig.KeyChans{3});
                            EEG = eeg_checkset( EEG );
                        else
                            disp('No .bdf file found');
                            return
                        end
                    case '.xdf'
                        FileToOpen = [Subject_Path  SUB{i} '.xdf'];
                        if isfile(FileToOpen)
                            EEG = pop_loadxdf(FileToOpen);
                            EEG = eeg_checkset( EEG );
                        else
                            disp('No .xdf file found');
                            return
                        end
                end
            end
            
        else
            % not the first file to open. So grab the last file that we
            % processed.
            FileToOpen = [Subject_Path  SUB{i} '.set'];
            if isfile(FileToOpen)
                % shorten FileToOpen to relative address for EEGlab commands
                FileToOpen = [SUB{i} DataConfig.LastSuffix{1}];
                EEG = pop_loadset( 'filename', FileToOpen, 'filepath', Subject_Path);
            else
                MessageForUser = ['No relevant ' DataConfig.LastSuffix{1} ' file found'];
                disp(MessageForUser );
                return
            end
        end
        
        % earlier call before more flexible opening command included.
        % FileToOpen = [Subject_Path  SUB{i} '.bdf'];
        
        % Take the .bdf file and make it a .set file (Oren added this)
        % Token reference is channel 40 for input, and then will reference
        % properly in a moment. Files must be named: "23.bdf" on the way in
        % and will be named "17.set" on the way out.
        
        % If file are e.g. "P17.bdf" then just adjust SUB cell array with
        % appropriate values.
        
        if DataConfig.DownSample{1} == EEG.srate
            % already at the desired sample rate so leave it be.
        else
            % Downsample from the recorded sampling rate e.g. 2048 Hz to e.g. 256 Hz
            % to speed data processing (automatically applies the appropriate
            % low-pass anti-aliasing filter)
            EEG = pop_resample(EEG, DataConfig.DownSample{1});
        end
        % save downsampled, or non-downsampled, as new data set. 
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname',[SUB{i} '_ds'],...
            'savenew',[Subject_Path SUB{i} '_ds.set'] ,'gui','off');
        
        % For PREP, only calculate the HEOG and VEOG (no mastoid/average reference reference)
        % choose relevant channel montage.
        switch DataConfig.RawFileType{1}
            case '.bdf'
                NoRefChanLbls = [Current_File_Path filesep 'SupportingDocs' ...
                    filesep 'ChannelsFor' num2str(DataConfig.TotalChannels{1}) '_NoRef_BDF.txt'];
            case '.xdf'
                NoRefChanLbls = [Current_File_Path filesep 'SupportingDocs' ...
                    filesep 'ChannelsFor' num2str(DataConfig.TotalChannels{1}) '_NoRef_XDF.txt'];
        end
        % apply that montage.
        EEG = pop_eegchanoperator( EEG, NoRefChanLbls);
        
        % Add channel location information corresponding to the 3-D coordinates of the electrodes based on 10-10 International System site locations
        EEG = pop_chanedit(EEG, 'lookup',[Current_File_Path filesep 'SupportingDocs' filesep DataConfig.ChanLocs{1}]);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3, 'setname', [SUB{i} '_ds_addChans'], ...
            'savenew', [Subject_Path SUB{i} '_ds_addChans.set'], 'gui', 'off');
        
        %% and now the PREP section
        % Used for basic cleaning, high pass filtering and detection and
        % rejection (and interpolation) of bad scalp channels.
        
        % 7 steps in total. Let's do them all.
        
        % notably, the parameters needed for each sub-process don't seem to
        % autopopulate with their default values. So let's give them all
        % their default values in advance here. If needed, the individual
        % execution lines below can overwritte as required.
        
        % This was when we were doing it piecemeal. Not now.
        PrepConfig = prepare_PREP_parameters(EEG);
        
        rawEEG = EEG; % save a copy of the dirty data for later comparison.
        
        % The massive all in one step, which I can't get to work.
        % [EEG, PREP_params, computationTimes] = prepPipeline(EEG, PrepConfig);
        
        % 2. Remove trend (high pass) temporarily to properly compute thresholds
        % we have some non-EEG channels, so better exclude them.
        % use default 0.5Hz cut-off (but written in case we need to change).
        detrendIn = PrepConfig.detrendIn;
        [EEG, detrendOut] = removeTrend(EEG, detrendIn);
        % EEG = pop_saveset( EEG, 'filename',[SUB{i} '_ds_detrend.set'],'filepath',Subject_Path);
        
        % 3. Remove line noise without committing to a filtering strategy
        % just the scalp channels, and multiples of 50Hz for Australia.
        % otherwise let them have it. Operates like a fancy notch filter.
        % Apparently fails without specifying sample frequency explicitly?
        lineNoiseIn = PrepConfig.lineNoiseIn;
        [EEG, lineNoiseOut] = cleanLineNoise(EEG, lineNoiseIn);
        % EEG = pop_saveset( EEG, 'filename',[SUB{i} '_ds_deNoise.set'],'filepath',Subject_Path);
        
        % 4. Robustly reference the signal relative to an estimate of the “true” average reference
        % AND simultaneously
        % 5. Detect and interpolate bad channels relative to this reference
        
        % We have non-EEG channels and the algorithm looks at the "other"
        % channels to determine what's noisy. So we should ban access to
        % other non-EEG channels or it'll perform poorly. As to what to
        % apply to the reference too, we should include EOGs and mastoids.
        referenceIn = PrepConfig.referenceIn;
        if DataConfig.PREP{1} == 1 % do the full PREP robust referencing.
            [EEG, referenceOut] = performReference(EEG, referenceIn);
        else % just use a regular average reference.
            AvgRef = mean(EEG.data(DataConfig.firstScalp:DataConfig.lastScalp,:),1);
            EEG.data = EEG.data - AvgRef;
        end
        % EEG = pop_saveset( EEG, 'filename',[SUB{i} '_ds_ref.set'],'filepath',Subject_Path);
        % 6. Produce reports if desired
        % can't get this to work for now as it's a product of the
        % PIPELINE and that doesn't seem to be working, so let's
        % just leave it and do the comparisons ourselves.
        
        %                 summary_filename = [DIR 'PREP_Summary.html'];
        %                 subject_filename = [Subject_Path filesep 'OtherData' filesep SUB{i} '_PREP_output.html'];
        %                 publishPrepReport(EEG, summary_filename, subject_filename, 1, true);
        % 7. Post process if desired (that's all the rest of our steps...)
        
        % Done with PREP.
        
        % no filtering performed by PREP pipeline so do it here.
        % ERP CORE default settings: (non-causal Butterworth impulse response function, 0.1 Hz half-amplitude cut-off, 12 dB/oct roll-off)
        
        % save immediately post PREP
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4, 'setname', [SUB{i} '_ds_addChans_PREP'], ...
            'savenew', [Subject_Path SUB{i} '_ds_addChans_PREP.set'], 'gui', 'off');
        
        % do simple FFT of pre-filtered data and save it for later drawing.
        tempvar = [];
        unfilt_FFT = zeros(DataConfig.TotalChannels{1}, length(abs(SimpleFFT(rawEEG.data(1,:)))));
        for ThisChan = DataConfig.firstScalp:DataConfig.lastScalp
            tempvar = abs(SimpleFFT(rawEEG.data(ThisChan,:)));
            unfilt_FFT(ThisChan,:) = tempvar;
        end
        unfilt_FFT = mean(unfilt_FFT,1);
        freqs = [0:(EEG.srate/2)/(length(unfilt_FFT)-1):EEG.srate/2];
        
        EEG  = pop_basicfilter( EEG,  DataConfig.KeyChans{1}:DataConfig.KeyChans{2} , 'Boundary', 'boundary', ...
            'Cutoff',  [DataConfig.HPfilter{1} DataConfig.LPfilter{1}], ...
            'Design', 'butter', 'Filter', 'bandpass', 'Order',  DataConfig.FiltOrder{1}, 'RemoveDC', 'on' );
        
        % save at the bandpass point.
                [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5, 'setname', [SUB{i} '_ds_addChans_PREP_bp'], ...
            'savenew', [Subject_Path SUB{i} '_ds_addChans_PREP_bp.set'], 'gui', 'off');

        % do simple FFT of post-filtered data and save it for later
        filt_FFT = zeros(DataConfig.TotalChannels{1}, length(abs(SimpleFFT(EEG.data(1,:)))));
        for ThisChan = DataConfig.firstScalp:DataConfig.lastScalp
            tempvar = abs(SimpleFFT(EEG.data(ThisChan,:)));
            filt_FFT(ThisChan,:) = tempvar;
        end
        filt_FFT = mean(filt_FFT,1);
        
        % if you need to do some debugging of drawing, save workspace now.
        % save('Debug_workspace.mat');
        
        if ~exist([Subject_Path filesep 'Figures' ])
            % create the relevant output subfolders
            mkdir([Subject_Path filesep 'Figures' ]); % Folder for figures
            mkdir([Subject_Path filesep 'OtherData' ]); % Folder for summary stats
        end
        
        % currently have an average reference. Lets add mastoid reference.
        % but first have to put Mastoids on the common reference if that
        % was used for scalp channels (when processed with PREP).
        if DataConfig.PREP{1} == 1
            CommonReference = referenceOut.referenceSignal;
            
            EEG.data(DataConfig.KeyChans{3},:) = EEG.data(DataConfig.KeyChans{3},:) - CommonReference;
            EEG.data(DataConfig.KeyChans{4},:) = EEG.data(DataConfig.KeyChans{4},:) - CommonReference;
            
            % and for a fair comparison we should put the raw data on the same
            % reference.
            for ThisChannel = 1:size(rawEEG.data,1)
                rawEEG.data(ThisChannel,:) = rawEEG.data(ThisChannel,:) - CommonReference;
            end
            
            % do a basic baselining of whole epoch solely for visualizing the
            % cleaning processes. Do it for raw and clean. Baseline = 30s.
            % raw
            if size(rawEEG.data,2) > 30*rawEEG.srate
                for k = 1:size(rawEEG.data,1)
                    rawEEG.data(k,:) = rawEEG.data(k,:) - mean(rawEEG.data(k,1:30*rawEEG.srate));
                end
            else % not enough data, so do instantaneous baseline.
                for k = 1:size(rawEEG.data,1)
                    rawEEG.data(k,:) = rawEEG.data(k,:) - rawEEG.data(k,1);
                end
            end
            
            % do a basic baselining of whole epoch solely for visualizing the
            % cleaning processes. Do it for raw and clean. Baseline = 30s.
            % clean
            if size(EEG.data,2) > 30*EEG.srate
                for k = 1:size(EEG.data,1)
                    EEG.data(k,:) = EEG.data(k,:) - mean(EEG.data(k,1:30*EEG.srate));
                end
            else % not enough data, so do instantaneous baseline.
                for k = 1:size(EEG.data,1)
                    EEG.data(k,:) = EEG.data(k,:) - EEG.data(k,1);
                end
            end
            
            % save the cleaning output (so we know which channels were
            % interpolated and why).
            CleaningFilename = [Subject_Path filesep 'OtherData' filesep SUB{i} '_PREPcleaned_output.mat' ];
            save(CleaningFilename , 'referenceOut');
        end
        
        % finally rereference to Mastoids if that's specified.
        if strcmp(DataConfig.ReReference{1}, 'Mastoid')
            MastoidRef_filename = ['ChannelsFor' num2str(DataConfig.TotalChannels{1}) '_addMastoidsPostPREP.txt'];
            EEG = pop_eegchanoperator( EEG, [Current_File_Path filesep 'SupportingDocs' filesep MastoidRef_filename]);
            % EEG = pop_saveset( EEG, 'filename',[SUB{i} '_ds_ref_bandpass_mastoids.set'],'filepath',Subject_Path);
        else % leave the optimizied common (average) reference in play.
        end
        
        % somehow doing this loses channel location information?
        % Add channel location information corresponding to the 3-D coordinates of the electrodes based on 10-10 International System site locations
        EEG = pop_chanedit(EEG, 'lookup',[Current_File_Path filesep 'SupportingDocs' filesep DataConfig.ChanLocs{1}]);
        % save the output.
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 6, 'setname', [SUB{i} '_ds_addChans_PREP_bp_refs'], ...
            'savenew', [Subject_Path SUB{i} '_ds_addChans_PREP_bp_refs.set'], 'gui', 'off');
        
        % old command doesn't udpate the name of the .set file internally.
        % EEG = pop_saveset( EEG, 'filename',[SUB{i} '_ds_PREP.set'],'filepath',Subject_Path);
        
        %% do the visualizations for each subject and save them.
        
        % do a quick plot of effect of bandpass filter first.
        figure;
        
        % prepare values for figures.
        lim_freqs = freqs(freqs<60);
        lim_freqs = lim_freqs(2:end); % exclude DC offset for scaling y-axis
        lim_filt_FFT = filt_FFT(freqs<60);
        lim_filt_FFT = lim_filt_FFT(2:end);
        lim_unfilt_FFT = unfilt_FFT(freqs<60);
        lim_unfilt_FFT = lim_unfilt_FFT(2:end);
        
        y_max = max([ quantile(lim_unfilt_FFT,0.95),  quantile(lim_filt_FFT,0.95)]);
        
        % plot the values.
        plot(lim_freqs, lim_unfilt_FFT,'-b',lim_freqs,lim_filt_FFT, '-r');
        % filtered plotted on top.
        if y_max > 0
            ylim([0 y_max]);
        else
            ylim([0 20]);
        end
        title('BandPass. Raw=B,Filt=R.')
        
        save2pdf([Subject_Path filesep 'Figures' filesep 'X1_' SUB{i} '_BPFilter.pdf']);
        
        if DataConfig.PREP{1} == 1
            % calculate time (in samples) for plotting
            times = [1:size(EEG.data, 2)];
            
            % divide time into four bins (with three boundaries)
            q1 = round(length(times)/4);
            q2 = round(length(times)/2);
            q3 = 3*q1;
            q4 = length(times);
            
            if ~isempty(referenceOut.badChannels.all)
                % there are some bad channels to report.
                GoodChannels = [ DataConfig.firstScalp : DataConfig.lastScalp ];
                BadChannels = referenceOut.badChannels.all(:);
                BadChannels = flip(BadChannels);
                
                % update the GoodChannels list to remove BadChannels
                for k = 1:length(BadChannels)
                    ThisBadChannel = BadChannels(k);
                    % remove each bad channel from list of good channels
                    GoodChannels(ThisBadChannel) = [];
                end
                
                % how many channels do we have now?
                NoOfBadChans = length(BadChannels);
                NoOfGoodChans = length(GoodChannels);
                
                % put all the good signals together
                GoodChannelData = EEG.data(GoodChannels, :);
                BadChannelData = referenceOut.badSignalsUninterpolated;
                % calculate an average signal
                GoodChannelData_mean = mean(GoodChannelData,1);
                BadChannelData_mean = mean(BadChannelData, 1);
                
                % open a bad channel figure to start plotting.
                % create a new figure.
                figure
                hold on
                for k = 1 : NoOfBadChans
                    plot(times, BadChannelData(k,:));
                end
                hold off
                title('BadChannelsRemoved')
                ylim([-500,500])
                save2pdf([Subject_Path filesep 'Figures' filesep 'X1_' SUB{i} '_BadChannels.pdf']);
                
                % create a figure of mean bad channels across time.
                figure
                plot(times, BadChannelData_mean);
                title('MeanOfBadChannels')
                ylim([-500,500])
                save2pdf([Subject_Path filesep 'Figures' filesep 'X1_' SUB{i} '_MeanOfBadChannels.pdf']);
                
            else % there were no bad channels to report.
                BadChannels = [];
                GoodChannels = [ DataConfig.firstScalp : DataConfig.lastScalp ];
                % how many channels do we have now?
                NoOfBadChans = length(BadChannels);
                NoOfGoodChans = length(GoodChannels);
                
                % put all the good signals together
                GoodChannelData = EEG.data(GoodChannels, :);
                BadChannelData = [];
                GoodChannelData_mean = mean(GoodChannelData,1);
                BadChannelData_mean = [];
                
                % no need to plot bad channel time series, so skip that.
                
                
                % no matter whether there were bad channels removed, still good to
                % plot the "clean data" metrics.
                % plot individual good signals uninterpolated
                figure %of number of chans interpolated.
                ChanReport = [NoOfGoodChans, NoOfBadChans];
                bar(ChanReport);
                xticklabels({'Good', 'Bad'});
                save2pdf([Subject_Path filesep 'Figures' filesep 'X1_' SUB{i} '_ChansRemoved.pdf']);
                close(gcf);
                
                figure
                hold on
                for k = 1 : NoOfGoodChans
                    plot(times,GoodChannelData(k,:));
                end
                hold off
                Baseline =  mean(mean(GoodChannelData));
                title('GoodChannelsRetained');
                ylim([Baseline-500,Baseline+500]);
                save2pdf([Subject_Path filesep 'Figures' filesep  'X1_' SUB{i} '_GoodChannels.pdf']);
                
                % create a figure of good channels across time.
                figure
                plot(times, GoodChannelData_mean);
                title('MeanOfGoodChannels')
                ylim([-500,500])
                save2pdf([Subject_Path filesep 'Figures' filesep 'X1_' SUB{i} '_MeanOfGoodChannels.pdf']);
                
                % no matter whether there were bad channels removed, still good to
                % plot the "clean data" metrics, including topographies.
                % and 4 topographies for clean, 4 for dirty.
                MeanDataPerChannel_raw_q1 = mean(rawEEG.data(DataConfig.firstScalp:DataConfig.lastScalp,1:q1),2);
                MeanDataPerChannel_clean_q1 = mean(EEG.data(DataConfig.firstScalp:DataConfig.lastScalp,1:q1),2);
                MeanDataPerChannel_raw_q2 = mean(rawEEG.data(DataConfig.firstScalp:DataConfig.lastScalp,q1+1:q2),2);
                MeanDataPerChannel_clean_q2 = mean(EEG.data(DataConfig.firstScalp:DataConfig.lastScalp,q1+1:q2),2);
                MeanDataPerChannel_raw_q3 = mean(rawEEG.data(DataConfig.firstScalp:DataConfig.lastScalp,q2+1:q3),2);
                MeanDataPerChannel_clean_q3 = mean(EEG.data(DataConfig.firstScalp:DataConfig.lastScalp,q2+1:q3),2);
                MeanDataPerChannel_raw_q4 = mean(rawEEG.data(DataConfig.firstScalp:DataConfig.lastScalp,q3+1:q4),2);
                MeanDataPerChannel_clean_q4 = mean(EEG.data(DataConfig.firstScalp:DataConfig.lastScalp,q3+1:q4),2);
                % import the chanlocs needed.
                ChanLocs = EEG.chanlocs(DataConfig.firstScalp:DataConfig.lastScalp);
                
                % plot those topographies
                PctThreshold = 95;
                MaxVal_clean = prctile(GoodChannelData_mean, PctThreshold); % 99th percentile of cleaned data set.
                MaxVal_raw = prctile(BadChannelData_mean, PctThreshold); % 99th percentile of raw data set.
                
                % use the largest scale, put everything on that.
                if MaxVal_raw > MaxVal_clean
                    caxis_both = [-1*MaxVal_raw MaxVal_raw];
                else % use max val from cleaned set.
                    caxis_both = [-1*MaxVal_clean MaxVal_clean];
                end
                % make sure goes low to high.
                caxis_both = sort(caxis_both);
                
                figure
                subplot(4,2,1)
                topoplot(MeanDataPerChannel_raw_q1, ChanLocs)
                colorbar
                caxis(caxis_both)
                title('Raw in Quartile1')
                
                subplot(4,2,2)
                topoplot(MeanDataPerChannel_raw_q2, ChanLocs)
                colorbar
                caxis(caxis_both)
                title('Raw in Quartile2')
                
                subplot(4,2,3)
                topoplot(MeanDataPerChannel_raw_q3, ChanLocs)
                colorbar
                caxis(caxis_both)
                title('Raw in Quartile3')
                
                subplot(4,2,4)
                topoplot(MeanDataPerChannel_raw_q4, ChanLocs)
                colorbar
                caxis(caxis_both)
                title('Raw in Quartile4')
                
                subplot(4,2,5)
                topoplot(MeanDataPerChannel_clean_q1, ChanLocs)
                colorbar
                caxis(caxis_both)
                title('Clean in Quartile1')
                
                subplot(4,2,6)
                topoplot(MeanDataPerChannel_clean_q2, ChanLocs)
                colorbar
                caxis(caxis_both)
                title('Clean in Quartile2')
                
                subplot(4,2,7)
                topoplot(MeanDataPerChannel_clean_q3, ChanLocs)
                colorbar
                caxis(caxis_both)
                title('Clean in Quartile3')
                
                subplot(4,2,8)
                topoplot(MeanDataPerChannel_clean_q4, ChanLocs)
                colorbar
                caxis(caxis_both)
                title('Clean in Quartile4')
                % and save the topographies all in one go.
                save2pdf([Subject_Path filesep 'Figures' filesep  'X1_' SUB{i} '_CleanedTopos.pdf']);
            end
        end
        close all % close off those figure windows and start again.
    end % of subject by subject loop
    
    % record the last process performed.
    DataConfig.LastProcess = cellstr('X1_PreProcess');
    DataConfig.LastSuffix = cellstr('_ds_addChans_PREP_bp_refs.set');

end