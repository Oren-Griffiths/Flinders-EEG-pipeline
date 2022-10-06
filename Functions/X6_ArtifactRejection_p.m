%Script #6
%Operates on individual subject data
%Uses the output from Script #5: Elist_Bin_Epoch.m
%This script loads the epoched EEG data file from Script #5, interpolates bad channels listed in Excel file Interpolate_Channels_MMN.xls, and performs artifact rejection to remove noisy segments
%of EEG and segments containing uncorrected residual eye movements using the parameters tailored to an individual subject's data listed in the corresponding Excel file for that artifact.

% additionally draws the average ERP pre- and post-AR for each person.

% and this is now the parallelized version.
function X6_ArtifactRejection_p(DataConfig,SUB, imageType)

% initialize
close all;
clearvars -except DataConfig SUB imageType;

    % Location of the main study directory
    DIR = fileparts(pwd);
    
    % location of preprocessing files.
    Current_File_Path = pwd;
    
    %Load the Excel file with the list of thresholds and parameters for identifying C.R.A.P. with the simple voltage threshold algorithm for each subject
    [ndata2, text2, alldata2] = xlsread([Current_File_Path filesep  'SupportingDocs' filesep 'AR_Parameters_for_SVT_CRAP']);
    
    %Load the Excel file with the list of thresholds and parameters for identifying eyeblinks during the stimulus presentation period (using the original non-ICA corrected VEOG signal) with the moving window peak-to-peak algorithm for each subject
    [ndata3, text3, alldata3] = xlsread([Current_File_Path filesep 'SupportingDocs' filesep 'AR_Parameters_for_MW_Blinks']);
    
    %*************************************************************************************************************************************
    
    %Loop through each subject listed in SUB
    for i = 1:length(SUB)
        
        %Open EEGLAB and ERPLAB Toolboxes
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        
        %Define subject path based on study directory and subject ID of current subject
        Subject_Path = [DIR filesep SUB{i} filesep];
        
        %Load the epoched EEG data file outputted from Script #5 in .set EEGLAB file format
        FileToOpen = [SUB{i} DataConfig.LastSuffix{1}];
        EEG = pop_loadset( 'filename', FileToOpen, 'filepath', Subject_Path);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [SUB{i} '_ds_PREP_ica_corr_cbip_elist_bins_epoch'], 'gui', 'off');
        
        %%
        %Identify segments of EEG with C.R.A.P. artifacts using the simple voltage threshold algorithm with the parameters in the Excel file for this subject
        DimensionsOfFile2 = size(alldata2);
        for j = 1:DimensionsOfFile2(1)
            if isequal(SUB{i},num2str(alldata2{j,1}))
                if isequal(alldata2{j,2}, 'default')
                    % Channels = ChannelsToTest;
                    Channels = [DataConfig.firstScalp:DataConfig.lastScalp];
                else
                    Channels = str2num(alldata2{j,2});
                end
                ThresholdMinimum = alldata2{j,3};
                ThresholdMaximum = alldata2{j,4};
                TimeWindowMinimum = alldata2{j,5};
                TimeWindowMaximum = alldata2{j,6};
            end
        end
        
        EEG  = pop_artextval( EEG , 'Channel',  Channels, 'Flag', [1 2], 'Threshold', [ThresholdMinimum ThresholdMaximum], 'Twindow', [TimeWindowMinimum  TimeWindowMaximum] );
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3, 'setname', [SUB{i} '_ds_PREP_ica_corr_cbip_elist_bins_epoch_ar'],...
             'savenew', [Subject_Path filesep SUB{i} '_ds_PREP_ica_corr_cbip_elist_bins_epoch_ar.set'], 'gui', 'off');

        %% apply blink filter if needed (and requested in config file).
        if DataConfig.RemoveBlinks{1} == 1
            DimensionsOfFile3 = size(alldata3);
            for j = 1:DimensionsOfFile3(1)
                if isequal(SUB{i},num2str(alldata3{j,1}))
                    Channel = alldata3{j,2};
                    Threshold = alldata3{j,3};
                    TimeWindowMinimum = alldata3{j,4};
                    TimeWindowMaximum = alldata3{j,5};
                    WindowSize = alldata3{j,6};
                    WindowStep = alldata3{j,7};
                end
            end
            
            EEG  = pop_artmwppth( EEG , 'Channel',  Channel, 'Flag', [1 4], 'Threshold', Threshold, 'Twindow', [TimeWindowMinimum  TimeWindowMaximum], 'Windowsize', WindowSize, 'Windowstep', WindowStep );
            pop_saveset(EEG, 'filename', [SUB{i} '_ds_PREP_ica_corr_cbip_elist_bins_epoch_ar'], 'filepath', Subject_Path); 
%             [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4, 'setname', [SUB{i} '_ds_PREP_ica_corr_cbip_elist_bins_epoch_ar'],...
%                    'savenew', [Subject_Path filesep SUB{i} '_ds_PREP_ica_corr_cbip_elist_bins_epoch_ar.set'], 'gui', 'off');
        end % checking whether to apply blink filter.
        
        %% plot the number of trials rejected, which channels, etc.
        % Oren wrote this bit.
        data = EEG.data(DataConfig.firstScalp:DataConfig.lastScalp, :,:);
        chanlocs = EEG.chanlocs(DataConfig.firstScalp:DataConfig.lastScalp);
        SVT = ThresholdMaximum;
        % Subject_Path = [fileparts(pwd) '\102\'];
        if strcmp(imageType,'none')
        else
            DrawARfigs(data,SVT,chanlocs, Subject_Path, imageType);
        end
       
    end   %End subject loop

end
%*************************************************************************************************************************************
