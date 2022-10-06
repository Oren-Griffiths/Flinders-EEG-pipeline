
function CreateIndividualARfiles(DataConfig)
% little bit of code to make the files that ERP CORE processing stream
% needs in order to do AR (including at the pre-ICA stage). This code
% generates a file with the same values for every subject. If you want to
% dial in the specific parameters per person, you'll need to turn this off
% and make the file manually (or just make the file by calling this and
% then adjust manually).

Current_File_Path = pwd;
outFilePath = [Current_File_Path filesep 'SupportingDocs'];

%% create the ICA prep values first.
NoOfSubs = length(DataConfig.SUB);

% initialize cell array of correct size
ICAprep = cell(NoOfSubs+1, 4);
% give it the right headers.
headers = {'SUBID', 'AmpthValue', 'WindowValue', 'StepValue'};
ICAprep(1,1:4) = headers;
% add in subject numbers
ICAprep(2:NoOfSubs+1,1) = DataConfig.SUB';
% find the right values and add them in.
AmpthValue = 500;
WindowValue = 2000;
StepValue = 50;
ICAvalues = {AmpthValue, WindowValue, StepValue};
% add in the values by row.
for ThisRow = 2:(NoOfSubs+1)
    ICAprep(ThisRow,2:4) = ICAvalues;
end
% output that cell into an excel file.
% make sure that any existing file is scrubbed.
tempFilename = [outFilePath filesep 'ICA_Prep_Values.xlsx'];
delete tempFilename; % won't crash if absent, will just give a command line warning.
writecell(ICAprep,tempFilename);

%% AR_Parameters_for_MW_Blinks.xlsx
% initialize cell array of correct size
AR_blinks = cell(NoOfSubs+1, 7);
% add in the appropriate headers.
headers = {'SUBID', 'Channel(s)', 'Threshold', 'Time Window Minimum', ...
    'Time Window Maximum', 'Window Size', 'Window Step'};
AR_blinks(1,1:7) = headers;
% add in subject numbers
AR_blinks(2:NoOfSubs+1,1) = DataConfig.SUB';
% find the right values and add them in.
VEOGchan = DataConfig.rawVEOG; % the VEOG channel specified in DataConfig.
threshold = 200; % probably sensible to have this determined in DataConfig too as it can change.
% check to see if the epoch extends 25ms before trigger onset.
if DataConfig.EpochMin{1} < -25
    timeWinMin = -25;
else
    timeWinMin = DataConfig.EpochMin{1};
end
% same, but checking that the epoch goes for at least 225ms.
if DataConfig.EpochMax{1} > 225
    timeWinMax = 225;
else
    timeWinMax = DataConfig.EpochMax{1};
end
winSize = 175;
winStep = 10;
AR_blinkValues = {VEOGchan, threshold, timeWinMin, timeWinMax, winSize, winStep };
% add in the values by row.
for ThisRow = 2:(NoOfSubs+1)
    AR_blinks(ThisRow,2:7) = AR_blinkValues;
end
% output that cell into an excel file.
tempFilename = [outFilePath filesep 'AR_Parameters_for_MW_Blinks.xlsx'];
delete tempFilename;
writecell(AR_blinks,tempFilename);

%% AR_Parameters_for_SVT_CRAP.xlsx
% initialize cell array of correct size
AR_SVT = cell(NoOfSubs+1, 6);
% add in the appropriate headers.
headers = {'SUBID', 'Channel(s)', 'Threshold Minimum', 'Threshold Maximum', ...
    'Time Window Minimum', 'Time Window Maximum'};
AR_SVT(1,1:6) = headers;
% add in subject numbers
AR_SVT(2:NoOfSubs+1,1) = DataConfig.SUB';
% find the right values and add them in.
% channels to assess is always 'default' (all channels) by default.
threshMin = -1*DataConfig.BadChanThreshold{1};
threshMax = DataConfig.BadChanThreshold{1};
timeWinMin = DataConfig.EpochMin{1}; % the whole epoch by default.
timeWinMax =  DataConfig.EpochMax{1};
AR_SVT_values = {'default', threshMin, threshMax, timeWinMin, timeWinMax};
% add in the values by row.
for ThisRow = 2:(NoOfSubs+1)
    AR_SVT(ThisRow,2:6) = AR_SVT_values;
end
% output that cell into an excel file.
tempFilename = [outFilePath filesep 'AR_Parameters_for_SVT_CRAP.xlsx'];
delete tempFilename;
writecell(AR_SVT,tempFilename);

%% finally, initialize an observed ICA component file too.
if DataConfig.ResetICAcomponents == 1
    % leave it blank except for row numbers per subject.
    ICAcomps =  cell(NoOfSubs+1, 2);
    % add in the appropriate headers.
    headers = {'SUBID', 'Components'};
    ICAcomps(1,1:2) = headers;
    % add in subject numbers
    ICAcomps(2:NoOfSubs+1,1) = DataConfig.SUB';
    % output that cell into an excel file.
    writecell(ICAcomps,[outFilePath filesep 'ICA_components.xlsx']);
else % check that there is actually a file to spare.
    if exist([outFilePath filesep 'ICA_components.xlsx']) > 0
        % if there's a file, just leave the ICA components alone.
    else % but if there's not, write a blank one anyway.
        % leave it blank except for row numbers per subject.
        ICAcomps =  cell(NoOfSubs+1, 2);
        % add in the appropriate headers.
        headers = {'SUBID', 'Components'};
        ICAcomps(1,1:2) = headers;
        % add in subject numbers
        ICAcomps(2:NoOfSubs+1,1) = DataConfig.SUB';
        % output that cell into an excel file.
        writecell(ICAcomps,[outFilePath filesep 'ICA_components.xlsx']);
    end
end

% create a participant mask file that defaults to including everyone.
% use the ICAcomps variable to do it because it's very similar..
ICAcomps(2:end,2) = {1};
% and initialize and write the variable.
tempFilename = [outFilePath filesep 'PIDmask.xlsx'];
delete tempFilename;
writecell(ICAcomps,tempFilename);

end