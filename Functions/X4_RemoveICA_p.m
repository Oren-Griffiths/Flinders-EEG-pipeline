%Script #4
%Operates on individual subject data
%Uses the output from Script #3: Run_ICA.m
%This script loads the outputted semi-continuous EEG data file containing the ICA weights from Script #3, loads the list of ICA component(s) from the ICA_Components_MMN.xlsx Excel file, and removes the component(s) from the EEG.
%Note that if ICA weights were re-computed on the data, the component(s) to remove will need to be updated in the Excel file to match the new components (see Script #3: Run_ICA.m for further details).

function X4_RemoveICA_p(DataConfig, SUB)

% initialize
close all;
clearvars -except DataConfig SUB;

    % Location of the main study directory
    DIR = fileparts(pwd)
    
    % location of preprocessing files.
    Current_File_Path = pwd;
    
    % initialize an output variable to track components removed.
    DataConfig.ComponentsRemoved = {};
    %**********************************************************************************************************************************************************************
    
    %Loop through each subject listed in SUB
    for i = 1:length(SUB)
        
        %Open EEGLAB and ERPLAB Toolboxes
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        
        %Define subject path based on study directory and subject ID of current subject
        Subject_Path = [DIR filesep SUB{i} filesep];
        
        %Load the continuous EEG data file containing the ICA weights outputted from Script #3 in .set EEGLAB file format
        FileToOpen = [SUB{i} DataConfig.LastSuffix{1}];
        EEG = pop_loadset( 'filename', FileToOpen, 'filepath', Subject_Path);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[SUB{i} '_ds_addChans_PREP_bp_refs_event_icaWeighted'], 'gui','off');
        
        % grab the adjusted VEOG (Fp1) data for output figure
        rawVEOG = EEG.data(DataConfig.Fp1,:);
        
        if isfield(EEG.etc, 'ic_classification')
            % this field will only exist if ICLabel has been applied. 
            % Parameters and workflow taken from Makato's blog.
            % https://sccn.ucsd.edu/wiki/Makoto's_preprocessing_pipeline#High-pass_filter_the_data_at_1-Hz_.28for_ICA.2C_ASR.2C_and_CleanLine.29.2805.2F18.2F2022_updated.29
            % Perform IC rejection using ICLabel scores, not dipole fitting.
            % Much more conservative criterion for removing components.
            % Instead of 70% likelihood of brain (Makato method), I've gone for
            % brain activity being the most likely classification (50%).
            % A more conservative option would be to just remove eyes. The
            % approach chosen needs to be specified.
            
            ic_threshold = 0.5; % 50% makes sense here, as must be most likely classification because
            % metric of ICLabel output is posterior probability (e.g. sum to 1).
            brainIdx  = find(EEG.etc.ic_classification.ICLabel.classifications(:,1) >= ic_threshold);
            eyeIdx = find(EEG.etc.ic_classification.ICLabel.classifications(:,3) >= ic_threshold);
            %
            % insert a default behaviour, and make consistent with config files
            % that don't use an ICA mode function.
            % default is remove eyes.
            if isfield(DataConfig, 'ICAmode')
                % would be stored as a cell array
                if iscell(DataConfig.ICAmode)
                    % a choice exists. Use that.
                    ICAmode = DataConfig.ICAmode{1};
                elseif isempty(DataConfig.ICAmode)
                    disp('DataConfig missing ICAmode property. Will just remove eye components.');
                    ICAmode = 'removeEyes';
                end
            else % earlier config files don't have that field. Just remove eyes.
                disp('DataConfig missing ICAmode property. Will just remove eye components.');
                ICAmode = 'removeEyes';
            end
            
            % save the good and bad channel indices.
            if strcmp(ICAmode, 'keepBrain')
                save([Subject_Path 'OtherData' filesep 'BrainComponentsKept.mat'], 'brainIdx');
                EEG = pop_subcomp(EEG, brainIdx, 0, 1); % updates IClabel fields too. 
            elseif strcmp(ICAmode, 'removeEyes')
                save([Subject_Path 'OtherData' filesep 'EyeComponentsRejected.mat'], 'eyeIdx');
                EEG = pop_subcomp(EEG, eyeIdx, 0, 0); % updates IClabel fields too. 
            else
                disp('Misspecified ICA behaviour. No components removed');
            end

        else % use the manual ICA correction process (adapted from ERPlab)

            %Load list of ICA component(s) corresponding to ocular artifacts from Excel file ICA_Components.xlsx
            [ndata, ~, alldata] = xlsread([Current_File_Path filesep 'SupportingDocs' filesep 'ICA_Components.xlsx']);
            MaxNumComponents = size(alldata, 2);
            for j = 1:size(alldata,1)
                if isequal(SUB{i}, num2str(alldata{j,1}))
                    NumComponents = 0;
                    for k = 2:MaxNumComponents
                        if ~isnan(alldata{j,k})
                            NumComponents = NumComponents+1;
                        end
                        Components = [alldata{j,(2:(NumComponents+1))}];
                    end
                end
            end
            
            % record which components were removed.
            DataConfig.ComponentsRemoved{i} = Components;
            save([Subject_Path 'OtherData' filesep 'ManualComponentsRej.mat'], 'Components');
            
            %Perform ocular correction by removing the ICA component(s) specified above
            EEG = pop_subcomp( EEG, [Components], 0, 0);
            
        end
        % % and save the IC_corrected data (if needed).
        % EEG = pop_saveset( EEG, 'filename', [SUB{i} '_ds_PREP_ica_corr.set'],'filepath', Subject_Path);

        %Create a bipolar HEOG channel (HEOG_left minus HEOG_right) and a bipolar VEOG channel (VEOG_lower minus FP2) from the ICA corrected data; the original uncorrected HEOG and VEOG channels are retained for later artifact detection procedures
                
        % add HEOG 
        switch DataConfig.TotalChannels{1} 
            case 32
                EEG.data(end+1,:) = EEG.data(36,:) - EEG.data(35,:);
            case 64
                EEG.data(end+1,:) = EEG.data(68,:) - EEG.data(67,:);
        end
        EEG.nbchan = size(EEG.data,1);
        if ~isempty(EEG.chanlocs)
            EEG.chanlocs(end+1).labels = 'corr_HEOG';
        end
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        eeglab redraw
        % add VEOG
                switch DataConfig.TotalChannels{1} 
            case 32
                 EEG.data(end+1,:) = EEG.data(37,:) - EEG.data(1,:);
            case 64
                EEG.data(end+1,:) = EEG.data(69,:) - EEG.data(1,:);
                end
        EEG.nbchan = size(EEG.data,1);
        if ~isempty(EEG.chanlocs)
            EEG.chanlocs(end+1).labels = 'corr_VEOG';
        end
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        eeglab redraw
        
        %Add channel location information corresponding to the 3-D coordinates of the electrodes based on 10-10 International System site locations
        % (or rather  re-add it, as we've changed the channels in this new file). 
        EEG = pop_chanedit(EEG, 'lookup',[Current_File_Path filesep 'SupportingDocs' filesep DataConfig.ChanLocs{1}]);
        EEG = pop_saveset(EEG, 'filename', [SUB{i} '_ds_PREP_ica_corr_cbip.set'],'filepath', Subject_Path);

        % and now grab the adjusted VEOG (Fp1) data.
        corrVEOG = EEG.data(DataConfig.Fp1,:);
        
        % draw and output the corrected VEOG against the raw VEOG
        time = [1:size(EEG.data, 2)];
        figure;
        
        plot(time, rawVEOG, time, corrVEOG);
        title('Fp1: RawVEOG = blue,CorrectedVEOG = orange');
        ylim([-1000, 1000]);
        % and save a picture.
        save2pdf([Subject_Path filesep 'Figures' filesep 'X4_' SUB{i} '_RawVsCorrected_VEOG.pdf']);

    end % End subject loop

end
%**********************************************************************************************************************************************************************
