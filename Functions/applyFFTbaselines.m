

% takes a single series of spectral powers (1 x N) and computes specral
% baseline (compares each point with N nearest neighbours, by calculating
% each sample as a Z score of its neighbours. 

% updated 17.09.21

function [signal_out] = applyFFTbaselines (varargin)

% usage is either:
% applyFFTbaseline (signal)
% applyFFTbaseline (signal, mode)
% applyFFTbaseline (signal, mode, halfwidth)
% applyFFTbaseline (signal, mode, halfwidth, excs)
% applyFFTbaseline (signal, mode, halfwidth, excs, ignoreNeighbour)

% signal = takes single vector (one series of spectral powers/amps)
% should be column vector. 

% mode is either 'none' or 'subtract' or 'Zscore' or 'SNR'
% 'none' = no subtraction.
% 'subtract' = subtract the mean of the baseline bins (linear subtraction)
% 'Zscore' = similar to 'subtract' but further normalized by SD of baseline
% 'SNR' = computes -log10(signal/noise) relative to baseline (units dB)
% default is 'SNR'.

% halfwidth = no of samples in either direction that forms the baseline.
% Defaults to 10 (so a total of 20 samples compared). 

% excs = no of samples in baseline range to exclude. Typically one
% maximum value is excluded (to minimize bias in baseline). Default is 1. 

% ignoreNeighbour: if == 1 then ignore immediatly adjacent bins.
% default is 1. 


% fill in missing input vals with defaults.
switch nargin
    case 5 
        signal = varargin{1};
        mode = varargin{2};
        halfwidth = varargin{3};
        excs = varargin{4};
        ignoreNeighbour = varargin{5};
    case 4 
        signal = varargin{1};
        mode = varargin{2};
        halfwidth = varargin{3};
        excs = varargin{4};
        ignoreNeighbour = 1;
    case 3
        signal = varargin{1};
        mode = varargin{2};
        halfwidth = varargin{3};
        excs = 1;
        ignoreNeighbour = 1;
    case 2 % default value for
        signal = varargin{1};
        mode = varargin{2};
        halfwidth = 10;
        excs = 1;
        ignoreNeighbour = 1;
    case 1
        signal = varargin{1};
        mode = 'SNR';
        halfwidth = 10;
        excs = 1;
        ignoreNeighbour = 1;
    case 0
        disp('Error: No signal entered into baseline function.');
        return
end

% be robust to orientation of signal vector
% makes any input a col vector.
if size(signal,2) > size (signal,1)
    signal = signal'; %
end

% make sure the requested mode is sensible
if strcmp(mode, 'none') || strcmp(mode,'subtract')|| strcmp(mode,'SNR') ||...
        strcmp(mode,'Zscore')
else
    disp('Error: Unexpected mode value entered');
    return
end


% figure out vector to work with and initialize output vector. 
NoOfSamples = length(signal);
signal_out = zeros(size(signal));

% if signal too short to do baseline correction, let user know to do
% better.
if NoOfSamples < (2*halfwidth)+1
    signal_out = [];
    disp('Too few samples to baseline correct');
    return
end

% 
if ignoreNeighbour == 1
    % need to start/end a little earlier or later (by one sample each)
    Startpoint = (halfwidth+2);
    Endpoint = (NoOfSamples-halfwidth-2);
else
    Startpoint = (halfwidth+1);
    Endpoint = (NoOfSamples-halfwidth-1);
end

for ThisSample = Startpoint:Endpoint
    
    % find halfwidth = N nearest neighbours for comparison
    if ignoreNeighbour == 1
        base_samples = [signal(ThisSample-(halfwidth+1):(ThisSample - 2)); ...
            signal((ThisSample+2):ThisSample+(halfwidth+1)) ];
    else
        % base_samples = signal(ThisSample-halfwidth:ThisSample+halfwidth);
        base_samples = [signal(ThisSample-(halfwidth):(ThisSample - 1)); ...
            signal((ThisSample+1):ThisSample+(halfwidth)) ];
    end
    
    % recursively remove the maximum N times (N = excs input parameter)
    for K = 1:excs
        [MaxVal, MaxIdx] = max(base_samples, [], 1);
        base_samples(MaxIdx) = [];
    end % of repeatedly removing max value.

    base_mean = mean(base_samples, 'omitnan');
    base_std = std(base_samples, 'omitnan');
    
    % calculate the mean of the remaining, unremoved vals.
    switch mode
        case 'none'
            % just report the mean.
            signal_out(ThisSample) = signal(ThisSample);
        case 'subtract'
            % subtract that mean (for standard subtraction)
            signal_out(ThisSample) = signal(ThisSample) - base_mean;
        case 'Zscore'
            % subtract that mean, and normalize (for Z scores).
            signal_out(ThisSample) = (signal(ThisSample) - base_mean)/base_std;
        case 'SNR'
            % calculate SNR relative to mean
            signal_out(ThisSample) = 10*log10(signal(ThisSample)/base_mean);
    end
    
    % collect and save the baseline mean value if it is first or last
    % sample (so can use that for the edge cases)
    if ThisSample == Startpoint
        early_base = base_mean;
        early_std = base_std;
    end
    
    if ThisSample == Endpoint
        late_base = base_mean;
        late_std = base_std;
    end
        
end

% and fix up those annoying edge cases.
% early edge
for ThisSample = 1:(Startpoint-1)
       % calculate the mean of the remaining, unremoved vals.
    switch mode
        case 'none'
            % just report the mean.
            signal_out(ThisSample) = signal(ThisSample);
        case 'subtract'
            % subtract that mean (for standard subtraction)
            signal_out(ThisSample) = signal(ThisSample) - early_base;
        case 'Zscore'
            % subtract that mean, and normalize (for Z scores).
            signal_out(ThisSample) = (signal(ThisSample) - early_base)/early_std;
        case 'SNR'
            % calculate SNR relative to mean
            signal_out(ThisSample) = 10*log10(signal(ThisSample)/early_base);
    end
end



% late edge
for ThisSample = Endpoint+1:NoOfSamples
        switch mode
        case 'none'
            % just report the mean.
            signal_out(ThisSample) = signal(ThisSample);
        case 'subtract'
            % subtract that mean (for standard subtraction)
            signal_out(ThisSample) = signal(ThisSample) - late_base;
        case 'Zscore'
            % subtract that mean, and normalize (for Z scores).
            signal_out(ThisSample) = (signal(ThisSample) - late_base)/late_std;
        case 'SNR'
            % calculate SNR relative to mean
            signal_out(ThisSample) = 10*log10(signal(ThisSample)/late_base);
    end
end

end