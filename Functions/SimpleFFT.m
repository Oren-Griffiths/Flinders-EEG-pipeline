

function OutVector = SimpleFFT(InVector)
% OG simpleFFT 15.07.20
    signal_FFT = fft(InVector) / length(InVector);
    FreqLength = floor(length(InVector)/2)+1;
    OutVector = horzcat(signal_FFT(1),2*abs(signal_FFT(2:FreqLength-1)),signal_FFT(FreqLength));
end