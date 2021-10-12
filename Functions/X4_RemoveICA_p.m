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
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[SUB{i} '_ds_PREP_ica_weighted'], 'gui','off');
        
        % grab the adjusted VEOG (Fp1) data for output figure
        rawVEOG = EEG.data(DataConfig.Fp1,:);
        
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
        
        %Perform ocular correction by removing the ICA component(s) specified above
        EEG = pop_subcomp( EEG, [Components], 0);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname',[SUB{i} '_ds_PREP_ica_corr'],'savenew', [Subject_Path SUB{i} '_ds_PREP_ica_corr.set'],'gui','off');
        
        %Create a bipolar HEOG channel (HEOG_left minus HEOG_right) and a bipolar VEOG channel (VEOG_lower minus FP2) from the ICA corrected data; the original uncorrected HEOG and VEOG channels are retained for later artifact detection procedures
        EEG = pop_eegchanoperator( EEG, [Current_File_Path filesep 'SupportingDocs' filesep DataConfig.AddCorrVEOG{1}]);
        
        %Add channel location information corresponding to the 3-D coordinates of the electrodes based on 10-10 International System site locations
        % (or rather  re-add it, as we've changed the channels in this new file). 
        EEG = pop_chanedit(EEG, 'lookup',[Current_File_Path filesep 'SupportingDocs' filesep DataConfig.ChanLocs{1}]);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3, 'setname', [SUB{i} '_ds_PREP_ica_corr_cbip'], 'savenew', [Subject_Path SUB{i} '_ds_PREP_ica_corr_cbip.set'], 'gui', 'off');
        
        % and now grab the adjusted VEOG (Fp1) data.
        corrVEOG = EEG.data(DataConfig.Fp1,:);
        
        % draw and output the corrected VEOG against the raw VEOG
        time = [1:size(EEG.data, 2)];
        figure
        
       plot(time, rawVEOG, time, corrVEOG);
       title('Fp1: RawVEOG = blue,CorrectedVEOG = orange');
       ylim([-1000, 1000]);
       % and save a picture.
       save2pdf([Subject_Path filesep 'Figures' filesep 'X4_' SUB{i} '_RawVsCorrected_VEOG.pdf']);

    end % End subject loop
    
%     DataConfig.LastProcess = cellstr('X4_RemoveICA');
%     DataConfig.LastSUB = SUB(i); % last participant processed.
%     DataConfig.LastSuffix = cellstr('_ds_PREP_ica_corr_cbip.set');

end
%**********************************************************************************************************************************************************************
