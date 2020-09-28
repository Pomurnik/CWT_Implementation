function  [X] = mystft(x,N,Fs,window,overlap)
    windowLength = window*Fs;
    j = 1;
    k = 1;
    while((k+windowLength-1)<N)
        n = k;
        tfragm = x(n:n+windowLength-1);
        hanning = hann(windowLength);
        xn = tfragm'.*hanning;
        FFT = abs(fft(xn));
        X(:,j) = FFT;
        j = j+1;
        k = k+windowLength-overlap*windowLength;
    end
end