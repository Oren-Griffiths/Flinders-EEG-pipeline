

% takes a single series of spectral powers (1 x N) and computes specral
% baseline (compares each point with N nearest neighbours, substract mean
% of those N nearest neighbours... then moves on).

function [signal_out] = applyFFTbaseline (varargin)

% usage is either:
% applyFFTbaseline (signal)
% applyFFTbaseline (signal, halfwidth)
% applyFFTbaseline (signal, halfwidth, excs)


% signal = takes single vector (one series of spectral powers/amps)
% halfwidth = no of samples in either direction taht forms the baseline.
% Defaults to 10 (so a total of 20 samples compared). 
% excs = no of samples in baseline range to exclude. Typically one
% maximum value is excluded. So default to 1. 

% fill in missing input vals with defaults.
switch nargin
    case 3
        signal = varargin{1};
        halfwidth = varargin{2};
        excs = varargin{3};
        % do nothing. fully specified.
    case 2 % default value for
        signal = varargin{1};
        halfwidth = varargin{2};
        excs = 1;
    case 1
        signal = varargin{1};
        halfwidth = 10;
        excs = 1;
    case 0
        disp('No signal entered into baseline function.');
        return
end

% be robust to orientation of signal vector
% makes any input a col vector.
if size(signal,2) > size (signal,1)
    signal = signal'; %
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
Startpoint = (halfwidth+1);
Endpoint = (NoOfSamples-halfwidth-1);

for ThisSample = Startpoint:Endpoint
    
    % find halfwidth = N nearest neighbours for comparison
    base_samples = signal(ThisSample-halfwidth:ThisSample+halfwidth);
    
    % recursively remove the maximum N times (N = excs input parameter)
    for K = 1:excs
        [MaxVal, MaxIdx] = max(base_samples, [], 1);
        base_samples(MaxIdx) = [];
    end % of repeatedly removing max value.
    
    % calculate the mean of the remaining, unremoved vals.
    base_mean = mean(base_samples);   
    % subtract that mean.
    signal_out(ThisSample) = signal(ThisSample) - base_mean;
    
    % collect and save the baseline mean value if it is first or last
    % sample (so can use that for the edge cases)
    if ThisSample == Startpoint
        early_base = base_mean;
    end
    
    if ThisSample == Endpoint
        late_base = base_mean;
    end
        
end

% and fix up those annoying edge cases.
% early edge
for ThisSample = 1:(Startpoint-1)
    signal_out(ThisSample) = signal(ThisSample)-early_base;
end
% late edge
for ThisSample = Endpoint+1:NoOfSamples
    signal_out(ThisSample) = signal(ThisSample)-late_base;
end

end