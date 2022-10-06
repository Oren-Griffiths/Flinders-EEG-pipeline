
% parameters

keyChans = 1:64; 
% for global field power

% keyChans = [11, 46, 47, 4, 39, 38, 37, 3, 36];
% for a set of 9 sensors centred on Fz
% {FC1, FC2, FCz, Fz, F1, F2, AFz, AF3, AF4};

keyHz = [4 7]; % theta is 4-7Hz.

timeWindow = [60 160]; % N1 peak is ~110ms (but wavelet window is ~500ms)

allConditions = {'B1(112)' , 'B2(114)' , 'B3(122)' , ...
    'B4(124)' , 'B5(132)' , 'B6(134)'};

% fixed values across all figures.
colour_max = 3;
colour_min = -1;
colour_max_itc = 0.35;
colour_min_itc = 0;

% global colour scheme
colScheme = 'parula';

%% header info in which we load in e.g. config information

% what's the relevant config file called?
ConfigFileName = 'WIMR_Config_TalkListenCued';

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
% shorten variable name
SUB = DataConfig.SUB;
load('chanlocs.mat'); % create variable "chanlocs"
load('cndSUBS.mat'); % loads in two variants of SUB variable, but with only 
% patients (SUB_patient) or only controls (SUB_control)
% and with noisy participant (pid=126) already removed. 

%% goal #1, plot  theta spatial topo during ERP peak (100-200ms) for patients
for thisPID = 1:length(SUB_patient)
    % loop through each needed file, load it, grab needed data, push it into
    % some holder variable.
    filename = ['TF_output' filesep SUB_patient{thisPID} '_TFdata.mat'];
    load(filename); % create variable tf_data
    
    for thisCND = 1:length(allConditions)
        % want 64 chans, average over timeWindow, average over keyHz
        time_idx = (tf_data.cond(thisCND).times >= timeWindow(1)) & ...
            (tf_data.cond(thisCND).times <= timeWindow(2)) ;
        freq_idx = (tf_data.cond(thisCND).freqs >= keyHz(1)) & ...
            (tf_data.cond(thisCND).freqs <= keyHz(2)) ;
        disp(['processing_cnd_' allConditions{thisCND} '_in_PID_' SUB_patient{thisPID} ])
        for thisChan = 1:64
            spatial_patient_ersp(thisPID, thisCND, thisChan) =  ...
                mean(tf_data.cond(thisCND).chan(thisChan).ersp(freq_idx,time_idx), 'all', 'omitnan');
            spatial_patient_itc(thisPID, thisCND, thisChan) =  ...
                mean(tf_data.cond(thisCND).chan(thisChan).itc(freq_idx,time_idx), 'all', 'omitnan');
        end % of channel by channel loop
    end % of condition by condition loop
    
end % of PID loop

% create output folders if needed 
if exist('TF_images', 'dir') == 7
else
    mkdir 'TF_images'
end

save('TF_images/spatial_patient_ersp.mat', 'spatial_patient_ersp');
save('TF_images/spatial_patient_itc.mat', 'spatial_patient_itc');

% plot each condition
% colour_max = max(spatial_patient_ersp(:));
% colour_min = min(spatial_patient_ersp(:));
% colour_max_itc = max(abs(spatial_patient_itc(:)));
% colour_min_itc = min(abs(spatial_patient_itc(:)));

for thisCND = 1:length(allConditions)
    % draw each figure for ersp
figure;
topoplot( squeeze(mean(spatial_patient_ersp(:,thisCND,1:64),1, 'omitnan')), chanlocs); 
colorbar;
caxis([colour_min, colour_max]);
colormap(colScheme);
% name and save each figure.
f = gcf;
f.Units = 'inches';
f.OuterPosition = [0.5 0.5 5.5 5.5]; % make the figure 5 inches in size.
fig_filename = ['TF_images' filesep 'Patient_ersp_' allConditions{thisCND} '_TargetHzTopoplot.png'];
disp(['Saving topoimage image ' fig_filename]);
exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
close(gcf);

% and output a spreadsheet of these values too. 
spreadsheetData = squeeze(spatial_patient_ersp(:,thisCND,:));
excelFilename = ['TF_images' filesep 'spatial_patient_ersp.xlsx'];
writematrix(spreadsheetData, excelFilename, 'Sheet', allConditions{thisCND});
writecell(SUB_patient', excelFilename, 'Sheet', 'SUBS');

% draw each figure for itc
figure;
topoplot( squeeze(mean(abs(spatial_patient_itc(:,thisCND,:)),1, 'omitnan')), chanlocs); 
colorbar;
colormap(colScheme);
caxis([colour_min_itc, colour_max_itc]);

% name and save each figure.
f = gcf;
f.Units = 'inches';
f.OuterPosition = [0.5 0.5 5.5 5.5]; % make the figure 5 inches in size.
fig_filename = ['TF_images' filesep 'Patient_itc_' allConditions{thisCND} '_TargetHzTopoplot.png'];
disp(['Saving topoimage image ' fig_filename]);
exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
close(gcf);

% and output a spreadsheet of these values too. 
spreadsheetData = squeeze(abs(spatial_patient_itc(:,thisCND,:)));
excelFilename = ['TF_images' filesep 'spatial_patient_itc.xlsx'];
writematrix(spreadsheetData, excelFilename, 'Sheet', allConditions{thisCND});
writecell(SUB_patient', excelFilename, 'Sheet', 'SUBS');
end


%% goal #2, plot  theta spatial topo during ERP peak (100-200ms) for patients
for thisPID = 1:length(SUB_ctrl)
    % loop through each needed file, load it, grab needed data, push it into
    % some holder variable.
    filename = ['TF_output' filesep SUB_ctrl{thisPID} '_TFdata.mat'];
    load(filename); % create variable tf_data
    
    for thisCND = 1:length(allConditions)
        % want 64 chans, average over timeWindow, average over keyHz
        time_idx = (tf_data.cond(thisCND).times >= timeWindow(1)) & ...
            (tf_data.cond(thisCND).times <= timeWindow(2)) ;
        freq_idx = (tf_data.cond(thisCND).freqs >= keyHz(1)) & ...
            (tf_data.cond(thisCND).freqs <= keyHz(2)) ;
        disp(['processing_cnd_' allConditions{thisCND} '_in_PID_' SUB_patient{thisPID} ])
        for thisChan = 1:64
            spatial_ctrl_ersp(thisPID, thisCND, thisChan) =  ...
                mean(tf_data.cond(thisCND).chan(thisChan).ersp(freq_idx,time_idx), 'all', 'omitnan');
            spatial_ctrl_itc(thisPID, thisCND, thisChan) =  ...
                mean(tf_data.cond(thisCND).chan(thisChan).itc(freq_idx,time_idx), 'all', 'omitnan');
        end % of channel by channel loop
    end % of condition by condition loop
    
end % of PID loop

% create output folders if needed 
if exist('TF_images', 'dir') == 7
else
    mkdir 'TF_images'
end

save('TF_images/spatial_ctrl_ersp.mat', 'spatial_ctrl_ersp');
save('TF_images/spatial_ctrl_itc.mat', 'spatial_ctrl_itc');

% plot each condition
% colour_max = max(spatial_ctrl_ersp(:));
% colour_min = min(spatial_ctrl_ersp(:));
% colour_max_itc = max(abs(spatial_ctrl_itc(:)));
% colour_min_itc = min(abs(spatial_ctrl_itc(:)));

for thisCND = 1:length(allConditions)
    % draw each figure for ersp
figure;
topoplot( squeeze(mean(spatial_ctrl_ersp(:,thisCND,1:64),1, 'omitnan')), chanlocs); 
colorbar;
colormap(colScheme);
caxis([colour_min, colour_max]);

% name and save each figure.
f = gcf;
f.Units = 'inches';
f.OuterPosition = [0.5 0.5 5.5 5.5]; % make the figure 5 inches in size.
fig_filename = ['TF_images' filesep 'ctrl_ersp_' allConditions{thisCND} '_TargetHzTopoplot.png'];
disp(['Saving topoimage image ' fig_filename]);
exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
close(gcf);

% and output a spreadsheet of these values too. 
spreadsheetData = squeeze(spatial_ctrl_ersp(:,thisCND,:));
excelFilename = ['TF_images' filesep 'spatial_ctrl_ersp.xlsx'];
writematrix(spreadsheetData, excelFilename, 'Sheet', allConditions{thisCND});
writecell(SUB_ctrl', excelFilename, 'Sheet', 'SUBS');

% draw each figure for itc
figure;
topoplot( squeeze(mean(abs(spatial_ctrl_itc(:,thisCND,:)),1, 'omitnan')), chanlocs); 
colorbar;
colormap(colScheme);
caxis([colour_min_itc, colour_max_itc]);

% name and save each figure.
f = gcf;
f.Units = 'inches';
f.OuterPosition = [0.5 0.5 5.5 5.5]; % make the figure 5 inches in size.
fig_filename = ['TF_images' filesep 'ctrl_itc_' allConditions{thisCND} '_TargetHzTopoplot.png'];
disp(['Saving topoimage image ' fig_filename]);
exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
close(gcf);

% and output a spreadsheet of these values too. 
spreadsheetData = squeeze(abs(spatial_ctrl_itc(:,thisCND,:)));
excelFilename = ['TF_images' filesep 'spatial_ctrl_itc.xlsx'];
writematrix(spreadsheetData, excelFilename, 'Sheet', allConditions{thisCND});
writecell(SUB_ctrl', excelFilename, 'Sheet', 'SUBS');
end


%% goal 3, grab mean time-freq plot for whole epoch, all freqs for patients
% first, patients
% keyChans = [11, 46, 47, 4, 39, 38, 37, 3, 36];

for thisPID = 1:length(SUB_patient)
    % loop through each needed file, load it, grab needed data, push it into
    % some holder variable.
    filename = ['TF_output' filesep SUB_patient{thisPID} '_TFdata.mat'];
    load(filename); % create variable tf_data
    
    for thisCND = 1:length(allConditions)
        % want just the key channel, all times, all freqs.
        disp(['processing_cnd_' allConditions{thisCND} '_in_PID_' SUB_patient{thisPID} ])
        for thisChan_idx = 1:length(keyChans)
            thisChan = keyChans(thisChan_idx);
            holder_ersp(thisChan_idx, :, :) = tf_data.cond(thisCND).chan(thisChan).ersp;
            holder_itc(thisChan_idx, :, :) = abs(tf_data.cond(thisCND).chan(thisChan).itc);
        end % of channel by channel loop
        % 4d: pid, cnd, time, freq.
        temporal_patient_ersp(thisPID, thisCND, :, :) =  ...
            mean(holder_ersp,1, 'omitnan');
        % 4d: pid, cnd, time, freq.
        temporal_patient_itc(thisPID, thisCND, :, :) =  ...
            mean(holder_itc,1, 'omitnan');
    end % of condition by condition loop
end % of PID loop

save('TF_images/temporal_patient_ersp.mat', 'temporal_patient_ersp');
save('TF_images/temporal_patient_itc.mat', 'temporal_patient_itc');

% plot each condition
% colour_max = max(temporal_patient_ersp(:));
% colour_min = min(temporal_patient_ersp(:));
% colour_max_itc = max(abs(temporal_patient_itc(:)));
% colour_min_itc = min(abs(temporal_patient_itc(:)));
times = tf_data.cond(1).times;
freqs = tf_data.cond(1).freqs;

for thisCND = 1:length(allConditions)
    % draw each figure for ersp
figure;
% 4d: pid, cnd, time, freq.
contourf(times,freqs, squeeze(mean(temporal_patient_ersp(:,thisCND,:,:), 1, 'omitnan'))  );
colorbar;
colormap(colScheme);
caxis([colour_min, colour_max]);

% font properties
ax = gca;
ax.FontSize = 16;
ax.FontWeight = 'bold';
ax.LineWidth = 3;

% name and save each figure.
f = gcf;
f.Units = 'inches';
f.OuterPosition = [0.5 0.5 5.5 5.5]; % make the figure 5 inches in size.
fig_filename = ['TF_images' filesep 'Patient_ersp_' allConditions{thisCND} '_timefreq.png'];
disp(['Saving timefreq image ' fig_filename]);
exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
close(gcf);

% and output a spreadsheet of these values too. 

% draw each figure for itc
figure;
contourf(times,freqs, squeeze(mean(temporal_patient_itc(:,thisCND,:,:), 1, 'omitnan'))  );
colorbar;
colormap(colScheme);
caxis([colour_min_itc, colour_max_itc]);

% font properties
ax = gca;
ax.FontSize = 16;
ax.FontWeight = 'bold';
ax.LineWidth = 3;

% name and save each figure.
f = gcf;
f.Units = 'inches';
f.OuterPosition = [0.5 0.5 5.5 5.5]; % make the figure 5 inches in size.
fig_filename = ['TF_images' filesep 'Patient_itc_' allConditions{thisCND} '_timefreq.png'];
disp(['Saving timefreq image ' fig_filename]);
exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
close(gcf);

% % and output a spreadsheet of these values too. 
% spreadsheetData = squeeze(abs(spatial_patient_itc(:,thisCND,:)));
% excelFilename = ['TF_images' filesep 'spatial_patient_itc.xlsx'];
% writematrix(spreadsheetData, excelFilename, 'Sheet', allConditions{thisCND});
end

%% goal 4, grab mean time-freq plot for whole epoch, all freqs for control

for thisPID = 1:length(SUB_ctrl)
    % loop through each needed file, load it, grab needed data, push it into
    % some holder variable.
    filename = ['TF_output' filesep SUB_ctrl{thisPID} '_TFdata.mat'];
    load(filename); % create variable tf_data
    
    for thisCND = 1:length(allConditions)
        % want just the key channel, all times, all freqs.
        disp(['processing_cnd_' allConditions{thisCND} '_in_PID_' SUB_ctrl{thisPID} ])
        for thisChan_idx = 1:length(keyChans)
            thisChan = keyChans(thisChan_idx);
            holder_ersp(thisChan_idx, :, :) = tf_data.cond(thisCND).chan(thisChan).ersp;
            holder_itc(thisChan_idx, :, :) = abs(tf_data.cond(thisCND).chan(thisChan).itc);
        end % of channel by channel loop
        % 4d: pid, cnd, time, freq.
        temporal_ctrl_ersp(thisPID, thisCND, :, :) =  ...
            mean(holder_ersp,1, 'omitnan');
        % 4d: pid, cnd, time, freq.
        temporal_ctrl_itc(thisPID, thisCND, :, :) =  ...
            mean(holder_itc,1, 'omitnan');
    end % of condition by condition loop
end % of PID loop

save('TF_images/temporal_ctrl_ersp.mat', 'temporal_ctrl_ersp');
save('TF_images/temporal_ctrl_itc.mat', 'temporal_ctrl_itc');

% plot each condition
% colour_max = max(temporal_ctrl_ersp(:));
% colour_min = min(temporal_ctrl_ersp(:));
% colour_max_itc = max(abs(temporal_ctrl_itc(:)));
% colour_min_itc = min(abs(temporal_ctrl_itc(:)));
times = tf_data.cond(1).times;
freqs = tf_data.cond(1).freqs;

for thisCND = 1:length(allConditions)
    % draw each figure for ersp
figure;
% 4d: pid, cnd, time, freq.
contourf(times,freqs, squeeze(mean(temporal_ctrl_ersp(:,thisCND,:,:), 1, 'omitnan'))  );
colorbar;
colormap(colScheme);
caxis([colour_min, colour_max]);

% font properties
ax = gca;
ax.FontSize = 16;
ax.FontWeight = 'bold';
ax.LineWidth = 3;

% name and save each figure.
f = gcf;
f.Units = 'inches';
f.OuterPosition = [0.5 0.5 5.5 5.5]; % make the figure 5 inches in size.
fig_filename = ['TF_images' filesep 'Ctrl_ersp_' allConditions{thisCND} '_timefreq.png'];
disp(['Saving timefreq image ' fig_filename]);
exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
close(gcf);

% and output a spreadsheet of these values too. 

% draw each figure for itc
figure;
contourf(times,freqs, squeeze(mean(temporal_ctrl_itc(:,thisCND,:,:), 1, 'omitnan'))  );
colorbar;
colormap(colScheme);
caxis([colour_min_itc, colour_max_itc]);

% font properties
ax = gca;
ax.FontSize = 16;
ax.FontWeight = 'bold';
ax.LineWidth = 3;

% name and save each figure.
f = gcf;
f.Units = 'inches';
f.OuterPosition = [0.5 0.5 5.5 5.5]; % make the figure 5 inches in size.
fig_filename = ['TF_images' filesep 'Ctrl_itc_' allConditions{thisCND} '_timefreq.png'];
disp(['Saving timefreq image ' fig_filename]);
exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
close(gcf);

% % and output a spreadsheet of these values too. 
% spreadsheetData = squeeze(abs(spatial_patient_itc(:,thisCND,:)));
% excelFilename = ['TF_images' filesep 'spatial_patient_itc.xlsx'];
% writematrix(spreadsheetData, excelFilename, 'Sheet', allConditions{thisCND});
end