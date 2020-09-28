function  [x,psiWTF, psiT] = mycwt(sig,Wn,t,T,fs)
    sigFFT = fft(sig);
    psiMat = zeros(1, length(sig));
        
    s_max = 100;
    for s=1:s_max  
        [PSI,t] = mywavelet(s,T,fs,Wn);
        psiWTF(:,s) = PSI;
        psiT(:,s) = t;
        psiFFT = fft(PSI);
        psiMat = [psiMat; psiFFT];
        
    end
    psiMat = rot90(psiMat(2:end,:),2);
    psiSize = size(psiMat);
    psiMat = [ psiMat(:,1:floor(psiSize(2)/2)), zeros(psiSize(1),ceil(psiSize(2)/2)) ];
        
    cFFT = sigFFT.*psiMat;
    c = ifft(cFFT,[],2);
    cSize = size(c);
    x = [c(:,floor(cSize(2)/2):end ), c(:,1:floor(cSize(2)/2))];
end