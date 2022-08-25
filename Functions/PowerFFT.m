

function OutVector = PowerFFT(InVector)
    % OG PowerFFT 22.08.22. FFT, but returns power.   
    signal_FFT = fft(InVector) / length(InVector);
    FreqLength = floor(length(InVector)/2)+1;
    OutVector = horzcat(abs(signal_FFT(1)).^2,abs(signal_FFT(2:FreqLength-1)).^2,abs(signal_FFT(FreqLength)).^2);
end