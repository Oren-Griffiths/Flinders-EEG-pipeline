

function X3b_PlotICAtopos_p(DataConfig, SUB)

% initialize
close all;
clearvars -except DataConfig SUB;


try
    % Location of the main study directory
    DIR = fileparts(pwd)
    
    % location of preprocessing files.
    Current_File_Path = pwd;
    
    
    %Loop through each subject listed in SUB
    for i = 1:length(SUB)
        
        %Open EEGLAB and ERPLAB Toolboxes
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        
        %Define subject path based on study directory and subject ID of current subject
        Subject_Path = [DIR filesep SUB{i} filesep];
        
        % record where we're up to in case of crash.
        DataConfig.CurrentSUB = SUB{i};
        
        % Load the semi-continuous EEG data file outputted from Script #2 in .set EEGLAB file format
        FileToOpen = [SUB{i} DataConfig.LastSuffix{1}];
        EEG = pop_loadset( 'filename', FileToOpen, 'filepath', Subject_Path);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [SUB{i} '_ds_PREP_ica_weighted'], 'gui', 'off');
        
        figure
        pop_topoplot(EEG, 0, [1:size(EEG.icaweights,1)],[SUB{i} '_ds_PREP_ica_weighted'], [round(sqrt(size(EEG.icaweights,1)))+1 round(sqrt(size(EEG.icaweights,1)))+1] ,0,'electrodes','on');
        save2pdf([Subject_Path filesep 'Figures' filesep 'X3b_' SUB{i} '_ICA_Weights.pdf']);
        close all
        
    end % End subject loop
    
    DataConfig.LastProcess = cellstr('X3b_PlotICAtopos');
    DataConfig.LastSUB = SUB(i); % last participant processed.
    DataConfig.LastSuffix = cellstr('_ds_PREP_ica_weighted.set');
    
catch ME
    display('Error in X3_RunICA. Workspace saved.');
    save('Debug_workspace.mat');
    rethrow(ME);
end % end of TRY loop
end