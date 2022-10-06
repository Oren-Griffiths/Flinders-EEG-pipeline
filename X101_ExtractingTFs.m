%% input parameters.

% which condition do we want to process?
% currently written to loop through all conditions.

% can be done with contains, so don't need all parts. 
% somtimes triger 55 can be recorded as "condition55"
allConditions = {'B1(' , 'B2(' , 'B3(' , ...
    'B4(' , 'B5(' , 'B6(', 'B7(', 'B8(', ...
    'B9(', 'B10(', 'B11(', 'B12('};

% which channels do you want to use? 
% currently written to loop through all scalp channels. 
keyChans = 1:32; 

%% header structure grabs file and config data

% what's the relevant config file called?
ConfigFileName = 'Config_Danielle_051022';

Current_File_Path = pwd;
addpath('Functions');
ConfigFilePath = [Current_File_Path filesep 'SupportingDocs' filesep ConfigFileName '.xlsx'];
Options = detectImportOptions(ConfigFilePath);

for k = 1:numel(Options.VariableTypes)
    Options.VariableTypes{k} = 'char';
end
DataConfig = table2struct(readtable(ConfigFilePath, Options));
DataConfig = adjustConfigData(DataConfig);

% and open eeglab to access the EEGlab functions
eeglab;
% just shorten variable name
SUB = DataConfig.SUB;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% manual override for troubleshooting.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:length(SUB)
    tic;
    PIDfolder = [fileparts(pwd) filesep SUB{k}];
    inputFile = [PIDfolder filesep SUB{k} '_ds_PREP_ica_corr_cbip_elist_bins_epoch_ar.set'];
    
    %% find each file and open it up
    EEG = pop_loadset('filename',  inputFile);
    
    % each file is missing channel information for some reason, so add it now.
    EEG = pop_chanedit(EEG, 'lookup',[Current_File_Path filesep 'SupportingDocs' filesep DataConfig.ChanLocs{1}]);
    
    %% separate  them out so each data set = 1 bin type.
    
    % first, select the epochs that we want. 
    % populate a list with of the time-locking event per epoch.
    
    for thisCond = 1:length(allConditions)
     cond2use = allConditions{thisCond};
        
    epochvect = cell(1,EEG.trials);
    keepTrials = zeros(1,EEG.trials);
    for i=1:EEG.trials
        [~,t0] = min(abs(cell2mat(EEG.epoch(i).eventlatency)));
        % epochvect(i) = EEG.epoch(i).eventtype{t0};
        epochvect{i} = EEG.epoch(i).eventtype{t0};
        % if strcmp(epochvect{i},cond2use)
        if contains(epochvect{i},cond2use)
            keepTrials(i) = 1;
        end
    end
    
    % limit our data to those epochs.
    data  = EEG.data(:,:,logical(keepTrials));

    
    %%  prepare parameters for repeated tf call.
    frames = EEG.pnts;
    tlimits = [EEG.times(1) EEG.times(end)];
    srate = EEG.srate;
    % cycles = [3 0.5]; % standard parameters for now.
    cycles = [1 0.5]; % allows us to get down into theta range.
    % cycles = 0; % FFT with Hanning tapering (no wavelets)
    
    % move these to header later.
    % keyChans = 1:64; 
    freqRange = [0 30];
    tempResolution = 200; % how many points to output on x-axis. (44s data, so 440 = 0.1s)
    baseline = [-1000 0]; % baseline period measured in epoch time (ms)
    
    % time-freq for data, per channel.
    % initialise some varialbes
    tf_data.cond(thisCond).times = [];
    tf_data.cond(thisCond).freqs = [];
    
    for thisChan = keyChans       
        
%         [ersp,itc,powbase,times,freqs,erspboot,itcboot,tfdata] = ...
%             newtimef(data(thisChan,:,:),EEG.pnts, [tlimits], 256, cycles, ...
%             'timesout', tempResolution,  'winsize', 128, ...
%             'plotersp', 'off', 'plotitc', 'off', 'trialbase', 'on');
        
        % changed window size back to default (as don't have the luxury)
        % but can buy some freq precision by increasing padratio. 
        
        [ersp,itc,powbase,times,freqs,erspboot,itcboot,tfdata] = ...
            newtimef(data(thisChan,:,:),EEG.pnts, [tlimits], 256, cycles, ...
            'timesout', tempResolution, 'padratio', 8, 'winsize', 80,  ...
            'plotersp', 'off', 'plotitc', 'off', 'trialbase', 'on');
        
        %     if you want to plot the outputs, then use code of this form.
        %     figure; contourf(times, freqs, ersp);
        %     figure; contourf(times, freqs, abs(itc));
        
        tf_data.PID = SUB{k};
        tf_data.cond(thisCond).chan(thisChan).lbl = EEG.chanlocs(thisChan).labels;
        tf_data.cond(thisCond).chan(thisChan).ersp = ersp;
        tf_data.cond(thisCond).chan(thisChan).itc = abs(itc);
        % check to see if you need to write before writing. 
        if isempty(tf_data.cond(thisCond).times) 
          tf_data.cond(thisCond).times = times;  
        end
        
        if isempty(tf_data.cond(thisCond).freqs) 
          tf_data.cond(thisCond).freqs = freqs;  
        end
        
        
    end % of channel by channel loop
end % of condition by condition loop

% output per person.
if exist('TF_output', 'dir') == 7
else
    mkdir 'TF_output'
end

outName = [pwd filesep 'TF_output' filesep SUB{k} '_TFdata.mat'];
save(outName, 'tf_data');
clear tf_data; % and then start over.

disp(['PID ' SUB{k} ' performed in ' num2str(toc)    ' seconds' ]);
end % of PID looping cycle

% and save the chanlocs structure for later
chanlocs = EEG.chanlocs;
save('chanlocs.mat', 'chanlocs');

