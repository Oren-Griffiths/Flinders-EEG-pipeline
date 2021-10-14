
% a little function that explicitly writes all the default values for the
% PREP functions. Not sure why, but the commands are not written to
% specify (or guess) defaults. So we'll populate a structure (of
% structures) here, called PrepConfig. Subfields are the structures needed
% for all of the individual commands.

function [PrepConfig, DataConfigOut] = prepare_PREP_parameters(EEG, DataConfigIn)

% now DataConfig is passed in, so no need for global.
% also requires DataConfig to be passed out as DataConfigOut at the end.
DataConfig = DataConfigIn; 

% detrendIn parameters
% this can be skipped because 'removeTrend' parses input variables nicely.
% But it is added for consistency.
% Make sure to detrend mastoids, even though they don't get interpolated
% or referenced
PrepConfig.detrendIn.detrendChannels = [DataConfig.firstScalp: DataConfig.lastScalp, ...
    DataConfig.KeyChans{3}, DataConfig.KeyChans{4}];
PrepConfig.detrendIn.detrendType = 'high pass';
PrepConfig.detrendIn.detrendCutoff = 1;
PrepConfig.detrendIn.detrendStepSize = 0.02;

% lineNoiseIn parameters
PrepConfig.lineNoiseIn.fPassBand = [0 , DataConfig.DownSample{1}];
PrepConfig.lineNoiseIn.Fs = DataConfig.DownSample{1};
PrepConfig.lineNoiseIn.fScanBandWidth = 2;
PrepConfig.lineNoiseIn.lineFrequencies = [50];
PrepConfig.lineNoiseIn.lineNoiseChannels = [DataConfig.firstScalp: DataConfig.lastScalp];
PrepConfig.lineNoiseIn.maximumIterations = 10;
PrepConfig.lineNoiseIn.p = 0.01;
PrepConfig.lineNoiseIn.pad = 0;
% some extra parameters that guide multitaper calculations.
PrepConfig.lineNoiseIn.taperBandWidth = 2;
PrepConfig.lineNoiseIn.taperWindowSize = 4;
PrepConfig.lineNoiseIn.taperWindowStep = 4;
PrepConfig.lineNoiseIn.tau = 100;

% referenceIn parameters
% evaluation channels must have a location in EEG.chanlocs, so best to
% avoid the external (EX) channels in this calculation.
PrepConfig.referenceIn.evaluationChannels = [DataConfig.firstScalp : DataConfig.lastScalp];
PrepConfig.referenceIn.rereference =  [DataConfig.firstScalp:DataConfig.lastScalp];
if strcmp(DataConfig.ReReference{1}, 'Mastoid')
    PrepConfig.referenceIn.referenceChannels = [DataConfig.KeyChans{3}, DataConfig.KeyChans{4}];
    PrepConfig.referenceIn.referenceType = 'specific'; % just does mastoids.
else % just do the regular rereferencing.
    % that is, a robust (interpolated) reference based on all scalp chans.
    PrepConfig.referenceIn.referenceType = 'robust';
    PrepConfig.referenceIn.referenceChannels = [DataConfig.firstScalp : DataConfig.lastScalp];
end
PrepConfig.referenceIn.interpolationOrder  = 'post-reference';% 'post-reference'; or 'post';
PrepConfig.referenceIn.meanEstimateType = 'median'; % 'median'
PrepConfig.referenceIn.channelLocations = EEG.chanlocs(DataConfig.firstScalp : DataConfig.lastScalp);
PrepConfig.referenceIn.channelInfo = EEG.chaninfo;
PrepConfig.referenceIn.srate = DataConfig.DownSample{1};
PrepConfig.referenceIn.samples = size(EEG.data,2);
PrepConfig.referenceIn.robustDeviationThreshold = 5;
PrepConfig.referenceIn.highfrequencyNoiseThreshold = 5;
PrepConfig.referenceIn.correlationWindowSeconds = 1;
% must have at least modest correlations between neighbouring channels.
% Could reduce for 32 channel montage, and will need to reduce or turn off
% for e.g. VR montage with far fewer channels. 
PrepConfig.referenceIn.correlationThreshold = 0.4;
PrepConfig.referenceIn.badTimeThreshold = 0.01;

% RANSAC is finding quite a few more channels than other criteria, and I
% understand it the least. So turning it off for now. 
PrepConfig.referenceIn.ransacOff = true(1); % true means NO_RANSAC.
PrepConfig.referenceIn.ransacSampleSize = 50;
PrepConfig.referenceIn.ransacChannelFraction = 0.25;
PrepConfig.referenceIn.ransacCorrelationThreshold = 0.75;
PrepConfig.referenceIn.ransacUnbrokenTime = 0.4;

PrepConfig.referenceIn.ransacWindowSeconds = 5;
PrepConfig.referenceIn.maxReferenceIterations = 4;
PrepConfig.referenceIn.reportingLevel = 'verbose';

DataConfigOut = DataConfig;

end