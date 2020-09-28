function  [PSI,t] = mywavelet(j,T,fs,Wn)
    s = (2^(1/20))^j; 
    t = s^2*(-T/2:1/fs:T/2);
   
    PSI = real(s^2*(pi^(-1/4)).*exp(1i.*Wn.*t).*exp(-(t.^2)./2));   
end