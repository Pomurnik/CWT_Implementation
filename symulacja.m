clear all
close all
clc

%% Parametry symulacji

Fs = 10000;
T = 10;
t = 0:1/Fs:T;
N = length(t);
M1 = [1 10.1];
M2 = [1 20000];
M3 = [1 900];
C1 = 10;
C2 = 1;
C3 = 10;
K1 = 200000;
K2 = 4000000;
K3 = 600000;

sim('model_3dof')

figure
subplot(3,1,1)
plot(tout,masa1)
title('Masa pierwsza')
xlabel('Czas [s]')
ylabel('Masa [kg]')
subplot(3,1,2)
plot(tout,masa2)
title('Masa druga')
xlabel('Czas [s]')
ylabel('Masa [kg]')
subplot(3,1,3)
plot(tout,masa3)
title('Masa trzecia')
xlabel('Czas [s]')
ylabel('Masa [kg]')

figure
subplot(3,1,1)
plot(tout,Res1)
title('Sygna³ odpowiedzi pierwszej masy')
xlabel('Czas [s]')
ylabel('Amplituda [-]')
subplot(3,1,2)
plot(tout,Res2)
title('Sygna³ odpowiedzi drugiej masy')
xlabel('Czas [s]')
ylabel('Amplituda [-]')
subplot(3,1,3)
plot(tout,Res3)
title('Sygna³ odpowiedzi trzeciej masy')
xlabel('Czas [s]')
ylabel('Amplituda [-]')

%% Fourier odpowiedzi na trzech masach

xf1 = abs(fft(Res1));
xf2 = abs(fft(Res2));
xf3 = abs(fft(Res3));

figure
subplot(3,1,1)
plot([0:N-1]/N*Fs, xf1)
xlim([0 500])
title('Transformata Fouriera sygna³u 1')
xlabel('Czêstotliwoœæ [Hz]')
ylabel('Amplituda')
subplot(3,1,2)
plot([0:N-1]/N*Fs, xf2)
xlim([0 500])
title('Transformata Fouriera sygna³u 2')
xlabel('Czêstotliwoœæ [Hz]')
ylabel('Amplituda')
subplot(3,1,3)
plot([0:N-1]/N*Fs, xf3)
xlim([0 500])
title('Transformata Fouriera sygna³u 3')
xlabel('Czêstotliwoœæ [Hz]')
ylabel('Amplituda')


%% Krótkookresowa transformata Fouriera

window = 2;
noverlap = 0.5;

mySpectrogram = mystft(Res1',N,Fs,window,noverlap);
mySpectrogram = mySpectrogram((1:length(mySpectrogram)/2),:);

figure
subplot(2,1,1);
imagesc([window/2 max(t)-window/2],[0 Fs/2],mySpectrogram)
colorbar
axis('xy')
title('Spektrogram w³asna implementacja')
xlabel('Czas [s]')
ylabel('Czêstotliwoœæ [Hz]')
ylim([0 100]);
subplot(2,1,2);
spectrogram(Res1',20000,10000,20000,Fs,'yaxis')
title('Spektrogram Matlab')
xlabel('Czas [s]')
ylabel('Czêstotliwoœæ [kHz]')
ylim([0 0.1]);

%% Transformata falkowa

Wn = 7;
[transformataF,psiWTF,psiT] = mycwt(Res1',Wn,t,T,Fs);

nrF = 1; 
figure
plot(psiT(:,nrF),psiWTF(:,nrF));
ylim([min(real(psiWTF(:,nrF))) max(real(psiWTF(:,nrF)))])
xlim([min(psiT(:,nrF)) max(psiT(:,nrF))])
if nrF == 1
    title("Falka matka");
    xlabel("Czas [s]");
    ylabel("Amplituda [-]");
elseif nrF > 1
    title("Falka córka");
    xlabel("Czas [s]");
    ylabel("Amplituda [-]");
end

figure
imagesc(t, [100 0], abs(transformataF));
axis('xy')
title('Transformata falkowa w³asna implementacja')
xlabel('Czas [s]')
ylim([0 100]);
figure
cwt(Res1',Fs)

